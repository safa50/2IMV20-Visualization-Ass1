/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package volvis;

import com.jogamp.opengl.GL;
import com.jogamp.opengl.GL2;
import com.jogamp.opengl.util.texture.Texture;
import com.jogamp.opengl.util.texture.awt.AWTTextureIO;
import gui.RaycastRendererPanel;
import gui.TransferFunction2DEditor;
import gui.TransferFunctionEditor;
import gui.TransferFunction2DEditor.TriangleWidget;
import java.awt.image.BufferedImage;
import util.TFChangeListener;
import util.VectorMath;
import volume.GradientVolume;
import volume.Volume;
import volume.VoxelGradient;

/**
 *
 * @author michel
 */
public class RaycastRenderer extends Renderer implements TFChangeListener {

    private Volume volume = null;
    private GradientVolume gradients = null;
    RaycastRendererPanel panel;
    TransferFunction tFunc;
    TransferFunctionEditor tfEditor;
    TransferFunction2DEditor tfEditor2D;

    public static final int SLICER_MODE = 0;
    public static final int MIP_MODE = 1;
    public static final int COMPOSITING_MODE = 2;
    public static final int TF2D_MODE = 3;

    private int mode = 0;
    private int step = 1;

    public boolean shading = false;

    // Set shading to true or false when shading checkbox is clicked
    public void setShading(boolean shading) {
        this.shading = shading;
    }

    // Set the rendering mode
    public void setMode(int mode) {
        this.mode = mode;
    }

    public RaycastRenderer() {
        panel = new RaycastRendererPanel(this);
        panel.setSpeedLabel("0");
    }

    public void setVolume(Volume vol) {
        System.out.println("Assigning volume");
        volume = vol;

        System.out.println("Computing gradients");
        gradients = new GradientVolume(vol);

        // set up image for storing the resulting rendering
        // the image width and height are equal to the length of the volume diagonal
        int imageSize = (int) Math.floor(Math.sqrt(vol.getDimX() * vol.getDimX() + vol.getDimY() * vol.getDimY()
                + vol.getDimZ() * vol.getDimZ()));
        if (imageSize % 2 != 0) {
            imageSize = imageSize + 1;
        }
        image = new BufferedImage(imageSize, imageSize, BufferedImage.TYPE_INT_ARGB);
        // create a standard TF where lowest intensity maps to black, the highest to white, and opacity increases
        // linearly from 0.0 to 1.0 over the intensity range
        tFunc = new TransferFunction(volume.getMinimum(), volume.getMaximum());

        // uncomment this to initialize the TF with good starting values for the orange dataset 
        //tFunc.setTestFunc();
        tFunc.addTFChangeListener(this);
        tfEditor = new TransferFunctionEditor(tFunc, volume.getHistogram());

        tfEditor2D = new TransferFunction2DEditor(volume, gradients);
        tfEditor2D.addTFChangeListener(this);

        System.out.println("Finished initialization of RaycastRenderer");
    }

    public RaycastRendererPanel getPanel() {
        return panel;
    }

    public TransferFunction2DEditor getTF2DPanel() {
        return tfEditor2D;
    }

    public TransferFunctionEditor getTFPanel() {
        return tfEditor;
    }

    short getVoxel(double[] coord) {

        if (coord[0] < 0 || coord[0] > volume.getDimX() || coord[1] < 0 || coord[1] > volume.getDimY()
                || coord[2] < 0 || coord[2] > volume.getDimZ()) {
            return 0;
        }

        int x = (int) Math.floor(coord[0]);
        int y = (int) Math.floor(coord[1]);
        int z = (int) Math.floor(coord[2]);

        return volume.getVoxel(x, y, z);
    }

    // Tri-linear interpolated voxel value
    short getTriVoxel(double[] coord) {

        if (coord[0] < 0 || coord[0] > volume.getDimX() || coord[1] < 0 || coord[1] > volume.getDimY()
                || coord[2] < 0 || coord[2] > volume.getDimZ()) {
            return 0;
        }

        double x = coord[0];
        double y = coord[1];
        double z = coord[2];

        int x0 = (int) Math.floor(x);
        int y0 = (int) Math.floor(y);
        int z0 = (int) Math.floor(z);

        int x1 = x0 + 1;
        int y1 = y0 + 1;
        int z1 = z0 + 1;

        if (x1 >= volume.getDimX() || y1 >= volume.getDimY() || z1 >= volume.getDimZ()) {
            return 0;
        }

        double alpha = x - x0;
        double beta = y - y0;
        double gamma = z - z0;

        short sx0 = volume.getVoxel(x0, y0, z0);
        short sx1 = volume.getVoxel(x1, y0, z0);
        short sx2 = volume.getVoxel(x0, y1, z0);
        short sx3 = volume.getVoxel(x1, y1, z0);
        short sx4 = volume.getVoxel(x0, y0, z1);
        short sx5 = volume.getVoxel(x1, y0, z1);
        short sx6 = volume.getVoxel(x0, y1, z1);
        short sx7 = volume.getVoxel(x1, y1, z1);

        short v = (short) ((1 - alpha) * (1 - beta) * (1 - gamma) * sx0
                + alpha * (1 - beta) * (1 - gamma) * sx1
                + (1 - alpha) * beta * (1 - gamma) * sx2
                + alpha * beta * (1 - gamma) * sx3
                + (1 - alpha) * (1 - beta) * gamma * sx4
                + alpha * (1 - beta) * gamma * sx5
                + (1 - alpha) * beta * gamma * sx6
                + alpha * beta * gamma * sx7);

        return v;
    }

    // Tri-linear interpolated gradient magnitude value
    double getTriGradientMag(double[] coord) {
        if (coord[0] < 0 || coord[0] > volume.getDimX() || coord[1] < 0 || coord[1] > volume.getDimY()
                || coord[2] < 0 || coord[2] > volume.getDimZ()) {
            return 0;
        }

        double x = coord[0];
        double y = coord[1];
        double z = coord[2];

        int x0 = (int) Math.floor(x);
        int y0 = (int) Math.floor(y);
        int z0 = (int) Math.floor(z);

        int x1 = x0 + 1;
        int y1 = y0 + 1;
        int z1 = z0 + 1;

        if (x1 >= volume.getDimX() || y1 >= volume.getDimY() || z1 >= volume.getDimZ()) {
            return 0;
        }

        double alpha = x - x0;
        double beta = y - y0;
        double gamma = z - z0;

        double sx0 = gradients.getGradient(x0, y0, z0).mag;
        double sx1 = gradients.getGradient(x1, y0, z0).mag;
        double sx2 = gradients.getGradient(x0, y1, z0).mag;
        double sx3 = gradients.getGradient(x1, y1, z0).mag;
        double sx4 = gradients.getGradient(x0, y0, z1).mag;
        double sx5 = gradients.getGradient(x1, y0, z1).mag;
        double sx6 = gradients.getGradient(x0, y1, z1).mag;
        double sx7 = gradients.getGradient(x1, y1, z1).mag;

        double v = ((1 - alpha) * (1 - beta) * (1 - gamma) * sx0
                + alpha * (1 - beta) * (1 - gamma) * sx1
                + (1 - alpha) * beta * (1 - gamma) * sx2
                + alpha * beta * (1 - gamma) * sx3
                + (1 - alpha) * (1 - beta) * gamma * sx4
                + alpha * (1 - beta) * gamma * sx5
                + (1 - alpha) * beta * gamma * sx6
                + alpha * beta * gamma * sx7);

        return v;
    }

    void slicer(double[] viewMatrix) {

        // clear image
        for (int j = 0; j < image.getHeight(); j++) {
            for (int i = 0; i < image.getWidth(); i++) {
                image.setRGB(i, j, 0);
            }
        }

        // vector uVec and vVec define a plane through the origin, 
        // perpendicular to the view vector viewVec
        double[] viewVec = new double[3];
        double[] uVec = new double[3];
        double[] vVec = new double[3];
        VectorMath.setVector(viewVec, viewMatrix[2], viewMatrix[6], viewMatrix[10]);
        VectorMath.setVector(uVec, viewMatrix[0], viewMatrix[4], viewMatrix[8]);
        VectorMath.setVector(vVec, viewMatrix[1], viewMatrix[5], viewMatrix[9]);

        // image is square
        int imageCenter = image.getWidth() / 2;

        double[] pixelCoord = new double[3];
        double[] volumeCenter = new double[3];
        VectorMath.setVector(volumeCenter, volume.getDimX() / 2, volume.getDimY() / 2, volume.getDimZ() / 2);

        // sample on a plane through the origin of the volume data
        double max = volume.getMaximum();
        TFColor voxelColor = new TFColor();

        for (int j = 0; j < image.getHeight(); j++) {
            for (int i = 0; i < image.getWidth(); i++) {
                pixelCoord[0] = uVec[0] * (i - imageCenter) + vVec[0] * (j - imageCenter)
                        + volumeCenter[0];
                pixelCoord[1] = uVec[1] * (i - imageCenter) + vVec[1] * (j - imageCenter)
                        + volumeCenter[1];
                pixelCoord[2] = uVec[2] * (i - imageCenter) + vVec[2] * (j - imageCenter)
                        + volumeCenter[2];

                int val = getVoxel(pixelCoord);

                // Map the intensity to a grey value by linear scaling
                voxelColor.r = val / max;
                voxelColor.g = voxelColor.r;
                voxelColor.b = voxelColor.r;
                voxelColor.a = val > 0 ? 1.0 : 0.0;  // this makes intensity 0 completely transparent and the rest opaque
                // Alternatively, apply the transfer function to obtain a color
                // voxelColor = tFunc.getColor(val);

                // BufferedImage expects a pixel color packed as ARGB in an int
                int c_alpha = voxelColor.a <= 1.0 ? (int) Math.floor(voxelColor.a * 255) : 255;
                int c_red = voxelColor.r <= 1.0 ? (int) Math.floor(voxelColor.r * 255) : 255;
                int c_green = voxelColor.g <= 1.0 ? (int) Math.floor(voxelColor.g * 255) : 255;
                int c_blue = voxelColor.b <= 1.0 ? (int) Math.floor(voxelColor.b * 255) : 255;
                int pixelColor = (c_alpha << 24) | (c_red << 16) | (c_green << 8) | c_blue;
                image.setRGB(i, j, pixelColor);
            }
        }

    }

    void mip(double[] viewMatrix) {
        // clear image
        for (int j = 0; j < image.getHeight(); j++) {
            for (int i = 0; i < image.getWidth(); i++) {
                image.setRGB(i, j, 0);
            }
        }

        // vector uVec and vVec define a plane through the origin, 
        // perpendicular to the view vector viewVec
        double[] viewVec = new double[3];
        double[] uVec = new double[3];
        double[] vVec = new double[3];
        VectorMath.setVector(viewVec, viewMatrix[2], viewMatrix[6], viewMatrix[10]);
        VectorMath.setVector(uVec, viewMatrix[0], viewMatrix[4], viewMatrix[8]);
        VectorMath.setVector(vVec, viewMatrix[1], viewMatrix[5], viewMatrix[9]);

        // image is square
        int imageCenter = image.getWidth() / 2;

        double[] pixelCoord = new double[3];
        double[] volumeCenter = new double[3];
        VectorMath.setVector(volumeCenter, volume.getDimX() / 2, volume.getDimY() / 2, volume.getDimZ() / 2);

        // sample on a plane through the origin of the volume data
        double max = volume.getMaximum();
        TFColor voxelColor = new TFColor();

        // When the user is toggeling with mouse, make it more interactive and faster
        // by increasing the step size.
        // If not toggeling, step size is one for good results
        if (this.interactiveMode) {
            this.step = 15;
        } else if (!this.interactiveMode) {
            this.step = 1;
        }

        for (int j = 0; j < image.getHeight(); j++) {
            for (int i = 0; i < image.getWidth(); i++) {

                int val = 0;
                for (int k = -volume.getDimZ(); k < volume.getDimZ(); k += step) {
                    pixelCoord[0] = uVec[0] * (i - imageCenter) + vVec[0] * (j - imageCenter)
                            + viewVec[0] * k + volumeCenter[0];
                    pixelCoord[1] = uVec[1] * (i - imageCenter) + vVec[1] * (j - imageCenter)
                            + viewVec[1] * k + volumeCenter[1];
                    pixelCoord[2] = uVec[2] * (i - imageCenter) + vVec[2] * (j - imageCenter)
                            + viewVec[2] * k + volumeCenter[2];

                    // Get the maximum voxel value
                    val = Math.max(val, getTriVoxel(pixelCoord));

                    // We use the Local Maximum Intensity Projection
                    if (val > 0.95 * max) {
                        break;
                    }
                }

                // Map the intensity to a grey value by linear scaling
                voxelColor.r = val / max;
                voxelColor.g = voxelColor.r;
                voxelColor.b = voxelColor.r;
                voxelColor.a = val > 0 ? 1.0 : 0.0;  // this makes intensity 0 completely transparent and the rest opaque
                // Alternatively, apply the transfer function to obtain a color
                // voxelColor = tFunc.getColor(val);

                // BufferedImage expects a pixel color packed as ARGB in an int
                int c_alpha = voxelColor.a <= 1.0 ? (int) Math.floor(voxelColor.a * 255) : 255;
                int c_red = voxelColor.r <= 1.0 ? (int) Math.floor(voxelColor.r * 255) : 255;
                int c_green = voxelColor.g <= 1.0 ? (int) Math.floor(voxelColor.g * 255) : 255;
                int c_blue = voxelColor.b <= 1.0 ? (int) Math.floor(voxelColor.b * 255) : 255;
                int pixelColor = (c_alpha << 24) | (c_red << 16) | (c_green << 8) | c_blue;
                image.setRGB(i, j, pixelColor);
            }
        }

    }

    void compositing(double[] viewMatrix) {
        // clear image
        for (int j = 0; j < image.getHeight(); j++) {
            for (int i = 0; i < image.getWidth(); i++) {
                image.setRGB(i, j, 0);
            }
        }

        // vector uVec and vVec define a plane through the origin, 
        // perpendicular to the view vector viewVec
        double[] viewVec = new double[3];
        double[] uVec = new double[3];
        double[] vVec = new double[3];
        VectorMath.setVector(viewVec, viewMatrix[2], viewMatrix[6], viewMatrix[10]);
        VectorMath.setVector(uVec, viewMatrix[0], viewMatrix[4], viewMatrix[8]);
        VectorMath.setVector(vVec, viewMatrix[1], viewMatrix[5], viewMatrix[9]);

        // image is square
        int imageCenter = image.getWidth() / 2;

        double[] pixelCoord = new double[3];
        double[] volumeCenter = new double[3];
        VectorMath.setVector(volumeCenter, volume.getDimX() / 2, volume.getDimY() / 2, volume.getDimZ() / 2);

        // sample on a plane through the origin of the volume data
        double max = volume.getMaximum();

        // When the user is toggeling with mouse, make it more interactive and faster
        // by increasing the step size.
        // If not toggeling, step size is one for good results
        if (this.interactiveMode) {
            this.step = 15;
        } else if (!this.interactiveMode) {
            this.step = 1;
        }

        for (int j = 0; j < image.getHeight(); j++) {
            for (int i = 0; i < image.getWidth(); i++) {

                TFColor previousColor = new TFColor();
                TFColor newColor = new TFColor();

                for (int k = -volume.getDimZ(); k < volume.getDimZ(); k += step) {
                    pixelCoord[0] = uVec[0] * (i - imageCenter) + vVec[0] * (j - imageCenter)
                            + viewVec[0] * k + volumeCenter[0];
                    pixelCoord[1] = uVec[1] * (i - imageCenter) + vVec[1] * (j - imageCenter)
                            + viewVec[1] * k + volumeCenter[1];
                    pixelCoord[2] = uVec[2] * (i - imageCenter) + vVec[2] * (j - imageCenter)
                            + viewVec[2] * k + volumeCenter[2];

                    int val = getTriVoxel(pixelCoord);
                    TFColor color = tFunc.getColor(val);

                    // Composite from back-to-front like in slides lecture 2
                    newColor.b = color.b * color.a + (1 - color.a) * previousColor.b;
                    newColor.g = color.g * color.a + (1 - color.a) * previousColor.g;
                    newColor.r = color.r * color.a + (1 - color.a) * previousColor.r;

                    newColor.a = (1 - color.a) * previousColor.a;

                    previousColor = newColor;

                }

                // BufferedImage expects a pixel color packed as ARGB in an int
                int c_alpha = 1 - newColor.a <= 1.0 ? (int) Math.floor((1 - newColor.a) * 255) : 255;
                int c_red = newColor.r <= 1.0 ? (int) Math.floor(newColor.r * 255) : 255;
                int c_green = newColor.g <= 1.0 ? (int) Math.floor(newColor.g * 255) : 255;
                int c_blue = newColor.b <= 1.0 ? (int) Math.floor(newColor.b * 255) : 255;
                int pixelColor = (c_alpha << 24) | (c_red << 16) | (c_green << 8) | c_blue;
                image.setRGB(i, j, pixelColor);
            }
        }
    }

    void tf2d(double[] viewMatrix) {
        // clear image
        for (int j = 0; j < image.getHeight(); j++) {
            for (int i = 0; i < image.getWidth(); i++) {
                image.setRGB(i, j, 0);
            }
        }

        // vector uVec and vVec define a plane through the origin, 
        // perpendicular to the view vector viewVec
        double[] viewVec = new double[3];
        double[] uVec = new double[3];
        double[] vVec = new double[3];
        VectorMath.setVector(viewVec, viewMatrix[2], viewMatrix[6], viewMatrix[10]);
        VectorMath.setVector(uVec, viewMatrix[0], viewMatrix[4], viewMatrix[8]);
        VectorMath.setVector(vVec, viewMatrix[1], viewMatrix[5], viewMatrix[9]);

        // image is square
        int imageCenter = image.getWidth() / 2;

        double[] pixelCoord = new double[3];
        double[] volumeCenter = new double[3];
        VectorMath.setVector(volumeCenter, volume.getDimX() / 2, volume.getDimY() / 2, volume.getDimZ() / 2);

        // sample on a plane through the origin of the volume data
        double max = volume.getMaximum();
        TFColor voxelColor = new TFColor();

        // When the user is toggeling with mouse, make it more interactive and faster
        // by increasing the step size.
        // If not toggeling, step size is one for good results
        if (this.interactiveMode) {
            this.step = 15;
        } else if (!this.interactiveMode) {
            this.step = 1;
        }

        //VoxelGradient gr = new VoxelGradient();        
        for (int j = 0; j < image.getHeight(); j++) {
            for (int i = 0; i < image.getWidth(); i++) {

                TFColor previousColor = new TFColor();
                TFColor newColor = new TFColor();

                for (int k = -volume.getDimZ(); k < volume.getDimZ(); k += step) {
                    pixelCoord[0] = uVec[0] * (i - imageCenter) + vVec[0] * (j - imageCenter)
                            + viewVec[0] * k + volumeCenter[0];
                    pixelCoord[1] = uVec[1] * (i - imageCenter) + vVec[1] * (j - imageCenter)
                            + viewVec[1] * k + volumeCenter[1];
                    pixelCoord[2] = uVec[2] * (i - imageCenter) + vVec[2] * (j - imageCenter)
                            + viewVec[2] * k + volumeCenter[2];

                    int fx = getTriVoxel(pixelCoord);
                    
                    //This for the out of bounds error
                    if (fx == 0) {
                        continue;
                    }

                    int fv = tfEditor2D.triangleWidget.baseIntensity;
                    double radius = tfEditor2D.triangleWidget.radius;
                    TFColor color = tfEditor2D.triangleWidget.color;

                    // These are the parameters for the Kniss method
                    double maxKnissGradMag = tfEditor2D.triangleWidget.maxKnissGradMag;
                    double minKnissGradMag = tfEditor2D.triangleWidget.minKnissGradMag;

                    VoxelGradient grad = gradients.getGradient((int) pixelCoord[0], (int) pixelCoord[1], (int) pixelCoord[2]);

                    //We use the Interpolated version to get the gradient magnitude.
                    //This is done because it gives smoother results than the standard: grad.mag.
                    double gradMag = getTriGradientMag(pixelCoord);

                    // Implement the Kniss constraint
                    if (gradMag > maxKnissGradMag || gradMag < minKnissGradMag) {
                        continue;
                    }

                    if (gradMag == 0 && fx == fv) {
                        voxelColor.a = color.a * 1.0;
                    } else if (gradMag > 0 && (fx - (radius * gradMag)) <= fv && fv <= (fx + (radius * gradMag))) {
                        voxelColor.a = color.a * (1.0 - ((1 / radius) * (Math.abs((fv - fx) / gradMag))));
                    } else {
                        voxelColor.a = 0.0;
                    }

                    // Phong Shading Model
                    if (shading) {
                        if (gradMag > 0.0 && voxelColor.a > 0.0) {
                            //Vector V = L = H, because the light is coming from us so flip view vector V.
                            double[] V = new double[]{-viewVec[0], -viewVec[1], -viewVec[2]};

                            // According to the paper of Levoy, vector L and H are normalized vectors.
                            // That's why we divide it by the length of V to normalize it.
                            double[] L = new double[3];
                            VectorMath.setVector(L, V[0] / VectorMath.length(V), V[1] / VectorMath.length(V), V[2] / VectorMath.length(V));

                            double[] H = new double[3];
                            VectorMath.setVector(H, V[0] / VectorMath.length(V), V[1] / VectorMath.length(V), V[2] / VectorMath.length(V));

                            double[] N = new double[3];
                            //Don't know why, but dividing by gradMag gives worse results for the shading
                            //so, divide by grad.mag.
                            VectorMath.setVector(N, grad.x / grad.mag, grad.y / grad.mag, grad.z / grad.mag);

                            //Parameters for Phong shading model
                            int alpha = 10;
                            double k_spec = 0.2;
                            double k_diff = 0.7;
                            double k_ambient = 0.1;

                            double LdotN = VectorMath.dotproduct(L, N);
                            double NdotH = VectorMath.dotproduct(N, H);

                            // Shading only works when the dotproducts are positive
                            if (LdotN > 0 && NdotH > 0) {
                                voxelColor.r = k_ambient + (color.r * k_diff * LdotN) + (k_spec * Math.pow(NdotH, alpha));
                                voxelColor.g = k_ambient + (color.g * k_diff * LdotN) + (k_spec * Math.pow(NdotH, alpha));
                                voxelColor.b = k_ambient + (color.b * k_diff * LdotN) + (k_spec * Math.pow(NdotH, alpha));
                            }
                        }

                        newColor.r = voxelColor.a * voxelColor.r + (1 - voxelColor.a) * previousColor.r;
                        newColor.g = voxelColor.a * voxelColor.g + (1 - voxelColor.a) * previousColor.g;
                        newColor.b = voxelColor.a * voxelColor.b + (1 - voxelColor.a) * previousColor.b;

                    } else {
                        newColor.r = voxelColor.a * color.r + (1 - voxelColor.a) * previousColor.r;
                        newColor.g = voxelColor.a * color.g + (1 - voxelColor.a) * previousColor.g;
                        newColor.b = voxelColor.a * color.b + (1 - voxelColor.a) * previousColor.b;
                    }

                    newColor.a = (1 - voxelColor.a) * previousColor.a;

                    previousColor = newColor;
                }

                // BufferedImage expects a pixel color packed as ARGB in an int
                int c_alpha = 1 - newColor.a <= 1.0 ? (int) Math.floor((1 - newColor.a) * 255) : 255;
                int c_red = newColor.r <= 1.0 ? (int) Math.floor(newColor.r * 255) : 255;
                int c_green = newColor.g <= 1.0 ? (int) Math.floor(newColor.g * 255) : 255;
                int c_blue = newColor.b <= 1.0 ? (int) Math.floor(newColor.b * 255) : 255;
                int pixelColor = (c_alpha << 24) | (c_red << 16) | (c_green << 8) | c_blue;
                image.setRGB(i, j, pixelColor);
            }
        }
    }

    private void drawBoundingBox(GL2 gl) {
        gl.glPushAttrib(GL2.GL_CURRENT_BIT);
        gl.glDisable(GL2.GL_LIGHTING);
        gl.glColor4d(1.0, 1.0, 1.0, 1.0);
        gl.glLineWidth(1.5f);
        gl.glEnable(GL.GL_LINE_SMOOTH);
        gl.glHint(GL.GL_LINE_SMOOTH_HINT, GL.GL_NICEST);
        gl.glEnable(GL.GL_BLEND);
        gl.glBlendFunc(GL.GL_SRC_ALPHA, GL.GL_ONE_MINUS_SRC_ALPHA);

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glDisable(GL.GL_LINE_SMOOTH);
        gl.glDisable(GL.GL_BLEND);
        gl.glEnable(GL2.GL_LIGHTING);
        gl.glPopAttrib();

    }

    @Override
    public void visualize(GL2 gl) {

        if (volume == null) {
            return;
        }

        drawBoundingBox(gl);

        gl.glGetDoublev(GL2.GL_MODELVIEW_MATRIX, viewMatrix, 0);

        long startTime = System.currentTimeMillis();

        // extended function to choose for different rendering methods
        switch (mode) {
            case SLICER_MODE:
                slicer(viewMatrix);
                break;
            case MIP_MODE:
                mip(viewMatrix);
                break;
            case COMPOSITING_MODE:
                compositing(viewMatrix);
                break;
            case TF2D_MODE:
                tf2d(viewMatrix);
                break;
            default:
                slicer(viewMatrix);
        }

        long endTime = System.currentTimeMillis();
        double runningTime = (endTime - startTime);
        panel.setSpeedLabel(Double.toString(runningTime));

        Texture texture = AWTTextureIO.newTexture(gl.getGLProfile(), image, false);

        gl.glPushAttrib(GL2.GL_LIGHTING_BIT);
        gl.glDisable(GL2.GL_LIGHTING);
        gl.glEnable(GL.GL_BLEND);
        gl.glBlendFunc(GL.GL_SRC_ALPHA, GL.GL_ONE_MINUS_SRC_ALPHA);

        // draw rendered image as a billboard texture
        texture.enable(gl);
        texture.bind(gl);
        double halfWidth = image.getWidth() / 2.0;
        gl.glPushMatrix();
        gl.glLoadIdentity();
        gl.glBegin(GL2.GL_QUADS);
        gl.glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
        gl.glTexCoord2d(0.0, 0.0);
        gl.glVertex3d(-halfWidth, -halfWidth, 0.0);
        gl.glTexCoord2d(0.0, 1.0);
        gl.glVertex3d(-halfWidth, halfWidth, 0.0);
        gl.glTexCoord2d(1.0, 1.0);
        gl.glVertex3d(halfWidth, halfWidth, 0.0);
        gl.glTexCoord2d(1.0, 0.0);
        gl.glVertex3d(halfWidth, -halfWidth, 0.0);
        gl.glEnd();
        texture.disable(gl);
        texture.destroy(gl);
        gl.glPopMatrix();

        gl.glPopAttrib();

        if (gl.glGetError() > 0) {
            System.out.println("some OpenGL error: " + gl.glGetError());
        }

    }
    private BufferedImage image;
    private double[] viewMatrix = new double[4 * 4];

    @Override
    public void changed() {
        for (int i = 0; i < listeners.size(); i++) {
            listeners.get(i).changed();
        }
    }
}
