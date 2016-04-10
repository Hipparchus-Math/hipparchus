package org.hipparchus.xdoc2md;

import java.io.IOException;
import java.nio.charset.Charset;
import java.nio.file.FileSystems;
import java.nio.file.FileVisitResult;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.SimpleFileVisitor;
import java.nio.file.attribute.BasicFileAttributes;
import java.util.List;


public class Converter {

    /**
     * Convert all files in arg[0] from xdoc to markdown and deposit the converted files into arg[1].
     *
     * @param args source and target paths
     */
    public static void main(String[] args) {
        final Path sourcePath = FileSystems.getDefault().getPath(args[0]);
        final Path targetPath = FileSystems.getDefault().getPath(args[1]);
        if (!sourcePath.toFile().isDirectory() || !targetPath.toFile().isDirectory() ) {
            System.out.println("Usage: Converter <source path> <target path>");
            System.exit(1);
        }
        try {
            convertDirectory(sourcePath, targetPath);
        } catch (Exception ex) {
            ex.printStackTrace();
            System.exit(1);
        }
    }

    /**
     * Converts xdoc file to markdown and places the converted document in target.
     *
     * @param xdoc source xdoc target
     * @param target destination markdown directory
     */
    private static void convertDocument(Path xdocPath, Path targetPath) {

        // Read the file into a List of lines
        Charset charset = Charset.defaultCharset();
        List<String> lines = null;
        try {
            lines = Files.readAllLines(xdocPath, charset);
        } catch (IOException ex) {
            ex.printStackTrace();
            System.exit(1);
        }

        //  Strip XML cruft
        lines = Utils.stripXML(lines);

        // Convert in-line <code> elements
        lines = Utils.convertCodeElements(lines);

        // Convert headings
        lines = Utils.convertHeadings(lines);

        // Convert links
        lines = Utils.convertLinks(lines);

        // Convert images
        lines = Utils.convertImageTags(lines);

        // Convert paragraphs
        lines = Utils.convertParagraphs(lines);

        // Convert tables
        lines = Utils.convertTables(lines);

        // Convert lists
        lines = Utils.convertLists(lines);

        // Convert source blocks and trim lines
        lines = Utils.convertSourceBlocksAndTrim(lines);

        // Write output file
        try {
            Files.write(targetPath, lines, charset);
        } catch (IOException ex) {
            ex.printStackTrace();
        }
    }

    /**
     * Converts all of the files in sourcePath to markdown and puts the
     * converted documents into targetPath. Currently does not traverse subdirectories.
     *
     * @param directory directory housing xdoc files.
     */
    private static void convertDirectory(final Path sourcePath, final Path targetPath) throws IOException {
        Files.walkFileTree(sourcePath, new SimpleFileVisitor<Path>() {
            @Override
            public FileVisitResult visitFile(Path sourceFile, BasicFileAttributes attrs)
                throws IOException {
                if (sourceFile.toFile().isDirectory())  {
                    throw new IllegalStateException("Directory to convert contains subdirectories - cannont handle this.");
                }
                if (!sourceFile.getFileName().toString().endsWith("xml")) {
                    return FileVisitResult.CONTINUE; // skip non-xml files
                }
                final Path destinationPath = FileSystems.getDefault().getPath(targetPath.toString(),
                                                   sourceFile.getFileName().toString().replace("xml", "md"));
                System.out.println("Converting " + sourceFile.toString() + " to " + destinationPath.toString());
                convertDocument(sourceFile, destinationPath);

                return FileVisitResult.CONTINUE;
            }
        });
    }
}
