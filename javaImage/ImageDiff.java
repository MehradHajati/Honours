import javax.imageio.ImageIO;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;

public class ImageDiff {

    public static void main(String[] args) {
        // File paths for the input images
        String imagePath1 = "image1.png";
        String imagePath2 = "image2.png";

        // Loading the images
        BufferedImage image1 = loadImage(imagePath1);
        BufferedImage image2 = loadImage(imagePath2);

        if (image1 != null && image2 != null) {
            // Compare images
            double difference = compareImages(image1, image2);
            System.out.println("Difference between the images: " + difference);
        } else {
            System.out.println("Failed to load images.");
        }
    }

    /**
     * Method to load in the images
     * @param imagePath Path of the png image
     * @return It will create a buffered image of the path passed into it
     */
    private static BufferedImage loadImage(String imagePath) {
        try {
            return ImageIO.read(new File(imagePath));
        } catch (IOException e) {
            e.printStackTrace();
            return null;
        }
    }

    /**
     * Method to compare the two images using the Euclidean distance for every correspoding pixel
     * @param image1 The first image
     * @param image2 The second image
     * @return Return the total difference
     */
    private static double compareImages(BufferedImage image1, BufferedImage image2) {
        // We can assume that the width and height of both of them will be the same
        int width = image1.getWidth();
        int height = image1.getHeight();

        double sumSquaredDiff = 0.0;

        // Iterate through pixels and calculate Euclidean distance
        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                int pixel1 = image1.getRGB(x, y) & 0xff; // Grayscale value of pixel in image1
                int pixel2 = image2.getRGB(x, y) & 0xff; // Grayscale value of pixel in image2
                // the bitwise operationg with 0xff is to discard the rgb value and only leave intensity value of the grayscale pixel

                // Calculate squared difference
                double diff = pixel1 - pixel2;
                sumSquaredDiff += diff * diff;
            }
        }

        // Calculate Euclidean distance
        return Math.sqrt(sumSquaredDiff);
    }
}