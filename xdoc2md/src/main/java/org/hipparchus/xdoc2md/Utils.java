package org.hipparchus.xdoc2md;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import org.jsoup.Jsoup;
import org.jsoup.nodes.Document;
import org.jsoup.nodes.Element;
import org.jsoup.select.Elements;

/**
 * Static methods used to convert xdoc fragments to markdown.
 */

public class Utils {

    final public static char QUOTE = Character.valueOf('"');
    final public static String START_SOURCE = "<source>";
    final public static String END_SOURCE = "</source>";
    final public static String START_CODE = "<code>";
    final public static String END_CODE = "</code>";
    final public static String START_BODY = "<body>";
    final public static String END_BODY = "</body>";
    final public static String END_DOCUMENT = "</document>";
    final public static String FOUR_SPACES = "    ";
    final public static String START_LINK = "<a href=";
    final public static String END_LINK = "</a>";
    final public static String SECTION_START = "<section name=";
    final public static String SUB_SECTION_START = "<subsection name=";
    final public static String SECTION_END = "</section>";
    final public static String SUB_SECTION_END = "</subsection>";
    final public static String START_P = "<p>";
    final public static String END_P = "</p>";
    final public static String TABLE_START = "<table";
    final public static String TABLE_END = "</table>";
    final public static String IMG_START = "<img src=";
    final public static String UL_START = "<ul>";
    final public static String UL_END = "</ul>";
    final public static String OL_START = "<ol>";
    final public static String OL_END = "</ol>";

    /**
     * Converts an html href link to markdown syntax
     * @param link converted link
     * @return converted link
     */
    public static String convertLink(String link)  {
        if (!link.startsWith("<a href=")) {
            throw new IllegalArgumentException("Malformed link: " + link);
        }
        char[] chars = link.toCharArray();
        String url = "";
        int pos = link.indexOf(QUOTE) + 1;
        while (chars[pos] != QUOTE) {
            url += chars[pos];
            pos++;
        }
        while (chars[pos] != '>') {
            pos++;
        }
        pos++;
        String label = "";
        while (chars[pos] != '<') {
            label += chars[pos];
            pos++;
        }
        return "[" + label + "]" + "(" + url + ")";
    }

    /**
     * Converts links in lines. Links may span multiple lines and there may be
     * multiple links per line.
     *
     * @param lines input lines
     * @return lines with links converted
     */
    public static List<String> convertLinks(List<String> lines) {
        ArrayList<String> out = new ArrayList<String>();
        int i = 0;
        while (i < lines.size()) {
            String newLine = lines.get(i);
            String outLine = "";
            if (newLine.contains(START_LINK)) {
                final int startLink = newLine.indexOf(START_LINK);
                if (newLine.contains(END_LINK) &&
                        newLine.lastIndexOf(END_LINK) > startLink) {  // link is on one line
                    String linkString =
                            newLine.substring(
                                              startLink,
                                              startLink + newLine.substring(startLink).indexOf(END_LINK) +
                                              END_LINK.length());
                    out.add(newLine.replace(linkString, convertLink(linkString)));
                    i++;
                    continue;
                }
                // multi-line link
                final int startLinkPos = newLine.indexOf(START_LINK);
                outLine = newLine.substring(0, startLinkPos);        // Text before the link
                String linkString = newLine.substring(startLinkPos); // Start of the link
                // Now need find the rest of the link in following line(s)
                int endPos = -1;
                int j = i + 1;
                String postLinkText = "";
                while (endPos < 0 && j < lines.size()) {
                    endPos = lines.get(j).indexOf(END_LINK);
                    if (endPos >= 0) {
                        linkString += lines.get(j).substring(0, endPos + END_LINK.length());
                        postLinkText = lines.get(j).substring(endPos + END_LINK.length());
                    } else {
                        linkString += lines.get(j);
                    }
                    j++;
                }
                i = j;
                out.add(outLine + convertLink(linkString));
                out.add(postLinkText);
            } else {
                out.add(newLine);
                i++;
            }
        }
        if (hasLinks(out)) {  // cheezy impl does not get multiple links per line on first pass...
            out = (ArrayList<String>) convertLinks(out);
        }
        return out;
    }

    /**
     * Checks to see if lines contains any links
     *
     * @param lines input lines to scan
     * @return true if there are any links in lines
     */
    public static boolean hasLinks(List<String> lines) {
        Iterator<String> iterator = lines.iterator();
        while (iterator.hasNext()) {
            if (iterator.next().contains(START_LINK)) {
                return true;
            }
        }
        return false;
    }

    /**
     * Converts img elements.  NOTE: cannot handle elements that span more than one line.
     *
     * @param lines input lines
     * @return lines with img elements converted to markdown format
     */
    public static List<String> convertImageTags(List<String> lines) {
        ArrayList<String> out = new ArrayList<String>();
        for (String newLine : lines) {
            String outLine = "";
            if (newLine.contains(IMG_START)) {
                final int startPos = newLine.indexOf(IMG_START);
                final int endPos = startPos + newLine.substring(startPos).indexOf(">");
                final String pre = newLine.substring(0, startPos);
                final String post = newLine.substring(endPos + 1);
                final String imgElt = newLine.substring(startPos, endPos + 1);
                final String[] values = imgElt.split("\"");
                String path = null;
                String altText = "";
                for (int i = 0; i < values.length; i++) {
                    final String curr = values[i].trim().toLowerCase();
                    if (curr.endsWith(IMG_START)) {
                        path = values[i+1];
                    }
                    if (curr.endsWith("alt=")) {
                        altText = values[i+1];
                    }
                }
                if (path == null) {
                    throw new IllegalArgumentException("Bad image element - cannot find path");
                }
                outLine = pre + "![" + altText + "](" + path + ")" + post;
            } else {
                outLine = newLine;
            }
            out.add(outLine);
        }
        return out;
    }

    /**
     * Converts <source> blocks to markdown code blocks.
     * Also trims leading and trailing spaces (other than the indents to indicate source blocks).
     *
     * @param lines input lines
     * @return output lines with code blocks converted
     */
    public static List<String> convertSourceBlocksAndTrim(List<String> lines) {
        boolean inBlock = false;
        Iterator<String> iterator = lines.iterator();
        ArrayList<String> out = new ArrayList<String>();
        while (iterator.hasNext()) {
            String nextLine = iterator.next();
            String outLine = "";
            if (nextLine.contains(START_SOURCE)) {
                nextLine = nextLine.trim();
                inBlock = true;
                if (nextLine.equals(START_SOURCE)) {
                    out.add(outLine);  // Add blank line to start block
                    continue;          // skip lines with just the start tag
                }
                if (nextLine.contains(END_SOURCE)) {  // one-line, inline code -> backticks
                    out.add(nextLine.replace(START_SOURCE, "`").replace(END_SOURCE, "`"));
                    inBlock = false;
                    continue;
                } else {
                    out.add(outLine);  // Add blank line to start block
                }
                // source tag is followed by code on the same line
                outLine = FOUR_SPACES + nextLine.substring(START_SOURCE.length());
            }
            if (inBlock) {
                if (nextLine.trim().equals(END_SOURCE)) {
                    inBlock = false;
                    continue;  // skip the line with just the end tag
                }
                if (outLine.equals("")) { // Not the first line in the block
                    outLine = FOUR_SPACES + nextLine;
                }
                if (outLine.contains(END_SOURCE)) {  // end tag at end of line
                    inBlock = false;
                    outLine = outLine.replace(END_SOURCE, "");
                }
            } else {
                outLine = nextLine.trim();  // not in source block, copy trimmed line
            }
            out.add(outLine);
        }
        return out;
    }

    /**
     * Converts <code> elements to markdown inline code (backticks)
     *
     * @param lines input lines
     * @return output lines with code elements converted
     */
    public static List<String> convertCodeElements(List<String> lines) {
        Iterator<String> iterator = lines.iterator();
        ArrayList<String> out = new ArrayList<String>();
        while (iterator.hasNext()) {
            String nextLine = iterator.next();
            out.add(nextLine.replace(START_CODE, "`").replace(END_CODE, "`"));
        }
        return out;
    }

    /**
     * Strips xml preamble up to and including <body> tag and end body, document tags
     *
     * @param lines input lines
     * @return lines between preamble and end body, doc tags
     */
    public static List<String> stripXML(List<String> lines) {
        Iterator<String> iterator = lines.iterator();
        ArrayList<String> out = new ArrayList<String>();
        boolean throatClearing = true;
        while (iterator.hasNext()) {
            final String nextLine = iterator.next();
            if (nextLine.contains(START_BODY)) {
                throatClearing = false;
                continue;
            }
            if (throatClearing) {
                continue;
            }
            final String trimmed = nextLine.trim();
            if (trimmed.equals(END_BODY) || trimmed.equals(END_DOCUMENT)) {
                continue;
            }
            out.add(nextLine);
        }
        return out;
    }

    /**
     * Converts section and subsection headings.  Converts subsection end tags to blank lines
     * and strips section end tag (no blank line at end).
     *
     * @param lines input lines
     * @return output lines with headings converted
     */
    public static List<String> convertHeadings(List<String> lines) {
        Iterator<String> iterator = lines.iterator();
        ArrayList<String> out = new ArrayList<String>();
        while (iterator.hasNext()) {
            String inLine = iterator.next();
            String outLine = null;
            if (inLine.contains(SECTION_START)) {
                inLine = inLine.trim();
                outLine = "# " + inLine.substring(inLine.indexOf(QUOTE) + 1, inLine.lastIndexOf(QUOTE));
            }
            if (inLine.contains(SUB_SECTION_START)) {
                inLine = inLine.trim();
                String[] chunks = inLine.split("\"");
                outLine = "## " + chunks[1];
            }
            if (inLine.contains(SUB_SECTION_END)) {
                outLine = "";  // blank line for subsection end
            }
            if (inLine.contains(SECTION_END)) {
                continue; // strip section end
            }
            if (outLine == null) {
                out.add(inLine);
            } else {
                out.add(outLine);
            }
        }
        return out;
    }

    /**
     * Converts paragraphs. Skips <p> lines and adds blank line in place of </p>.
     * Does not change indentation.
     *
     * @param lines input lines
     * @return lines with <p>'s converted
     */
    public static List<String> convertParagraphs(List<String> lines) {
        Iterator<String> iterator = lines.iterator();
        ArrayList<String> out = new ArrayList<String>();
        while (iterator.hasNext()) {
            String inLine = iterator.next();
            int lineType = 0;
            /*
             * 0 no tags
             * 1 <p>
             * 2 <p>blah
             * 3 <p>blah</p>
             * 4 </p>
             * 5 blah</p>
             */
            if (inLine.contains(START_P)) {
                inLine = inLine.trim();
                if (inLine.equals(START_P)) {
                    lineType = 1;
                }
                if (!inLine.startsWith(START_P)) {
                    throw new IllegalArgumentException("Cannot handle embedded <p>");
                }
                if (inLine.contains(END_P)) {
                    lineType = 3;
                }
                if (lineType == 0) {
                    lineType = 2;
                }
            } else {
                if (inLine.contains(END_P)) {
                    if (inLine.trim().equals(END_P)) {
                        lineType = 4;
                    } else {
                        lineType = 5;
                    }
                }
            }
            switch (lineType) {
                case 0:
                    out.add(inLine);
                    break;
                case 1:
                    // Skip <p>
                    break;
                case 2:
                    out.add(inLine.substring(START_P.length()));
                    break;
                case 3:
                    out.add(inLine.replace(START_P, "").replace(END_P, ""));
                    out.add("");
                    break;
                case 4:
                    out.add("");
                    break;
                case 5:
                    out.add(inLine.replace(END_P, ""));
                    out.add("");
                    break;
                default:
                    throw new IllegalStateException("Bad line type");
            }
        }
        return out;
    }

    /**
     * Converts tables to (extended) markdown syntax
     *
     * @param lines input lines
     * @return lines with tables converted
     */
    public static List<String> convertTables(List<String> lines) {
        final Iterator<String> iterator = lines.iterator();
        final ArrayList<String> out = new ArrayList<String>();
        String tableHTML = null;
        boolean inTable = false;
        while (iterator.hasNext()) {
            String inLine = iterator.next();
            if (inLine.contains(TABLE_START)) {
                tableHTML = "<html><body>" + inLine;
                inTable = true;
            } else if (inTable) {
                tableHTML += inLine;
            }
            if (inLine.contains(TABLE_END)) {
                if (!inLine.contains(TABLE_START)) {
                    tableHTML += inLine;
                }
                tableHTML += "</body></html>";

                // Should now have a valid html document containing the table to convert
                final Document doc = Jsoup.parse(tableHTML);
                Element table = doc.select("table").get(0); //select the first table.
                Elements rows = table.select("tr");

                // Create (mandatory) header line
                Elements headers = rows.get(0).select("th");
                String headerLine = "|";
                String dashLine = "";
                final boolean hasHeaders = !headers.isEmpty();
                if (!hasHeaders) {  // hacky, hacky - force first row to be th
                    rows.get(0).html(rows.get(0).html().replace("td>", "th>").replace("<td", "<th"));
                    headers = rows.get(0).select("th");
                }
                for (Element header : headers) {
                    headerLine += " " + header.html();
                    headerLine += " |";
                    dashLine += "| --- ";
                }
                out.add(headerLine);
                out.add(dashLine + "|");

                // Process rows
                for (Element row : rows) {
                    if (!row.select("th").isEmpty()) {
                        continue;  // skip header - already processed
                    }
                    Elements cells = row.select("td");
                    String rowLine = "|";
                    for (Element cell : cells) {
                        rowLine += " " + cell.html() + " |";
                    }
                    out.add(rowLine);
                    inTable = false;
                }
            } else {
                if (!inTable) {
                    out.add(inLine);
                }
            }
        }
        return out;
    }

    /**
     * Converts lists to markdown syntax
     *
     * @param lines input lines
     * @return lines with lists converted
     */
    public static List<String> convertLists(List<String> lines) {
        final Iterator<String> iterator = lines.iterator();
        final ArrayList<String> out = new ArrayList<String>();
        String listHtml = null;
        int ol = 0;
        int ul = 0;
        boolean inList = false;
        while (iterator.hasNext()) {
            final String inLine = iterator.next();
            if (inLine.contains(OL_START)) {
                inList = true;
                ol++;
                listHtml += inLine;
            } else if (inLine.contains(UL_START)) {
                inList = true;
                ul++;
                listHtml += inLine;
            } else if (inList) {
                listHtml += inLine;
            }
            if (inList && listHtml == null) {
                listHtml = "<html><body>";
            }
            if (inLine.contains(UL_END)) {
                ul--;
            }
            if (inLine.contains(OL_END)) {
                ol--;
            }
            if (!inList) {
                listHtml += inLine;
                out.add(inLine);
            }
            if (inList && ol + ul == 0) {
                listHtml += "</body></html>";
                final Document doc = Jsoup.parse(listHtml);
                final Element listElement = getContainedList(doc);
                if (listElement == null) {
                    throw new IllegalStateException("List not found.");
                }
                processList(listElement, out, 0);
                inList = false;
            }
        }
        return out;
    }

    /**
     * Process a JSOUP list Element, adding markdown to out.
     *
     * @param list List element
     * @param out  array of lines to add markdown to
     * @param level nesting level of the list
     */
    private static void processList(Element list, List<String> out, int level) {
        final int indentLevel = level * 4;
        boolean ol = list.tagName().equals("ol");
        if (!ol && !list.tagName().equals("ul")) {
            throw new IllegalStateException("Unrecognized list type " + list.tagName());
        }
        Elements listItems = list.children();
        for (Element item : listItems) {
            if (!item.tagName().equals("li")) {
                throw new IllegalStateException("Expecting <li> but got " + list.tagName());
            }
            final Element nestedList = getContainedList(item);
            if (nestedList != null) {
                processList(nestedList, out, level + 1);
            } else {
                out.add(spaces(indentLevel) + (ol ? "1. " : "* ") + item.html().trim());
            }
        }
        // Add blank line to end the list if the list is top-level
        if (level == 0) {
            out.add(" ");
        }
    }

    /**
     * Returns true if element is a list.
     *
     * @param element Element to examine
     * @return true iff Element is a list
     */
    private static boolean isList(Element element) {
        final String tagName = element.tagName();
        return tagName.equals("ol") || tagName.equals("ul");
    }

    /**
     * If element contains a list, return the first list it contains;
     * otherwise return null.
     *
     * @param element Element to examine
     * @return first embedded list if there is one; else null
     */
    private static Element getContainedList(Element element) {
        final Elements elements = element.getAllElements();
        Element listElement = null;
        for (int i = 1; i < elements.size(); i++) {
            if (isList(elements.get(i))) {
                listElement = elements.get(i);
                break;
            }
        }
        return listElement;
    }

    /**
     * Create a string of spaces of the given length.
     *
     * @param length length of returned string
     * @return a string of length spaces
     */
    private static String spaces(int length) {
        if (length == 0) {
            return "";
        }
        StringBuffer outputBuffer = new StringBuffer(length);
        for (int i = 0; i < length; i++){
            outputBuffer.append(" ");
        }
        return outputBuffer.toString();
    }

}

