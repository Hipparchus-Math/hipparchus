/*
 * Licensed to the Apache Software Foundation (ASF) under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The ASF licenses this file to You under the Apache License, Version 2.0
 * (the "License"); you may not use this file except in compliance with
 * the License.  You may obtain a copy of the License at
 *
 *      https://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

/*
 * This is not the original file distributed by the Apache Software Foundation
 * It has been modified by the Hipparchus project
 */
package org.hipparchus.samples;

import java.awt.Component;
import java.awt.Graphics2D;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.InputEvent;
import java.awt.event.KeyEvent;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;

import javax.imageio.ImageIO;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.KeyStroke;
import javax.swing.SwingUtilities;

/** Graphics utilities for examples.
 */
public class ExampleUtils {

    /** Empty constructor.
     * <p>
     * This constructor is not strictly necessary, but it prevents spurious
     * javadoc warnings with JDK 18 and later.
     * </p>
     * @since 3.0
     */
    public ExampleUtils() { // NOPMD - unnecessary constructor added intentionally to make javadoc happy
        // nothing to do
    }

    /** Display frame. */
    @SuppressWarnings("serial")
    public static class ExampleFrame extends JFrame {

        /** Empty constructor.
         * <p>
         * This constructor is not strictly necessary, but it prevents spurious
         * javadoc warnings with JDK 18 and later.
         * </p>
         * @since 3.0
         */
        public ExampleFrame() {
            // nothing to do
        }

        /**
         * Returns the main panel which should be printed by the screenshot action.
         * <p>
         * By default, it returns the content pane of this frame, but can be overriden
         * in case the frame has a global scroll pane which would cut off any offscreen content.
         *
         * @return the main panel to print
         */
        public Component getMainPanel() {
            return getContentPane();
        }
    }

    /** Display example.
     * @param frame frame to display
     */
    public static void showExampleFrame(final ExampleFrame frame) {
        Runnable r = new Runnable() {
            public void run() {
                JMenuItem screenshot = new JMenuItem("Screenshot (png)");
                screenshot.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_0, InputEvent.CTRL_DOWN_MASK));
                screenshot.addActionListener(new ActionListener() {
                    public void actionPerformed(ActionEvent ae) {
                        JFileChooser fileChooser = new JFileChooser(System.getProperty("user.dir"));
                        if (fileChooser.showSaveDialog(frame) == JFileChooser.APPROVE_OPTION) {
                          File file = fileChooser.getSelectedFile();
                          BufferedImage img = getScreenShot(frame.getMainPanel());
                          try {
                              // write the image as a PNG
                              ImageIO.write(img, "png", file);
                          } catch (IOException e) {
                              e.printStackTrace();
                          }
                        }
                    }
                });

                JMenuItem exit = new JMenuItem("Exit");
                exit.addActionListener(new ActionListener() {
                    public void actionPerformed(ActionEvent e) {
                        System.exit(0);
                    }
                });

                JMenu menu = new JMenu("File");
                menu.add(screenshot);
                menu.add(exit);
                JMenuBar mb = new JMenuBar();
                mb.add(menu);
                frame.setJMenuBar(mb);

                frame.setLocationRelativeTo(null);
                frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
                frame.setVisible(true);
            }
        };
        SwingUtilities.invokeLater(r);
    }

    private static BufferedImage getScreenShot(Component component) {
        BufferedImage image = new BufferedImage(component.getWidth(), component.getHeight(), BufferedImage.TYPE_INT_RGB);
        // call the Component's paint method, using the Graphics object of the image.
        component.paint(image.getGraphics());
        return image;
    }

    /** Resize an image.
     * @param originalImage original image
     * @param width desired width
     * @param height desired height
     * @param type type of the create image
     * @return resized image
     */
    public static BufferedImage resizeImage(BufferedImage originalImage, int width, int height, int type) {
        BufferedImage resizedImage = new BufferedImage(width, height, type);
        Graphics2D g = resizedImage.createGraphics();
        g.drawImage(originalImage, 0, 0, width, height, null);
        g.dispose();
        return resizedImage;
    }

}
