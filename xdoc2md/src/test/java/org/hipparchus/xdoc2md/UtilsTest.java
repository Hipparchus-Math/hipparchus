package org.hipparchus.xdoc2md;

import java.util.ArrayList;
import java.util.List;

import org.junit.Assert;
import org.junit.Test;

public class UtilsTest {

    @Test
    public void testConvertLink() {
        final String original = "<a href=" + Utils.QUOTE +
                "../apidocs/org.hipparchus/analysis/solvers/BisectionSolver.html" +
                Utils.QUOTE + ">Bisection</a>";
        final String converted = "[Bisection](../apidocs/org.hipparchus/analysis/solvers/BisectionSolver.html)";
        Assert.assertEquals(converted, Utils.convertLink(original));
    }

    @Test
    public void testConvertImg() {
        final List<String> lines = new ArrayList<String>();
        lines.add("    Deep and important commentary");
        lines.add("<img src=" + Utils.QUOTE + "../images/userguide/real_distribution_examples.png" +
                Utils.QUOTE + " alt=" + Utils.QUOTE + "Overview of continuous distributions" + Utils.QUOTE + "/>");
        lines.add(" and the deep stuff continues");
        final List<String> observed = Utils.convertImageTags(lines);
        Assert.assertEquals(lines.size(), observed.size());
        Assert.assertEquals(lines.get(0), observed.get(0));
        Assert.assertEquals(lines.get(2), observed.get(2));
        Assert.assertEquals("![Overview of continuous distributions](../images/userguide/real_distribution_examples.png)",
                            observed.get(1));
    }

    @Test
    public void testConvertImgNoAltText() {
        final List<String> lines = new ArrayList<String>();
        lines.add("<img src=" + Utils.QUOTE +
                  "../images/userguide/real_distribution_examples.png" + Utils.QUOTE + "/>");
        final List<String> observed = Utils.convertImageTags(lines);
        Assert.assertEquals(lines.size(), observed.size());
        Assert.assertEquals("![](../images/userguide/real_distribution_examples.png)",
                            observed.get(0));
    }

    @Test
    public void testConvertImgInline() {
        final List<String> lines = new ArrayList<String>();
        lines.add("  Boo <img src=" + Utils.QUOTE +
                  "../images/userguide/real_distribution_examples.png" + Utils.QUOTE + "/> hoo");
        final List<String> observed = Utils.convertImageTags(lines);
        Assert.assertEquals(lines.size(), observed.size());
        Assert.assertEquals("  Boo ![](../images/userguide/real_distribution_examples.png) hoo",
                            observed.get(0));
    }

    @Test
    public void testConvertCodeBlockTagOnFirstLine() {
        final List<String> inline = new ArrayList<String>();
        inline.add("    Deep and important commentary");
        inline.add("<source>UnivariateFunction function = // some user defined function object");
        inline.add("final double relativeAccuracy = 1.0e-12;");
        inline.add("final double absoluteAccuracy = 1.0e-8");
        inline.add("final int    maxOrder         = 5");
        inline.add("</source>");
        inline.add("And the saga continues  ");
        final List<String> expected = new ArrayList<String>();
        expected.add("Deep and important commentary");
        expected.add("");
        expected.add("    UnivariateFunction function = // some user defined function object");
        expected.add("    final double relativeAccuracy = 1.0e-12;");
        expected.add("    final double absoluteAccuracy = 1.0e-8");
        expected.add("    final int    maxOrder         = 5");
        expected.add("And the saga continues");
        final List<String> out = Utils.convertSourceBlocksAndTrim(inline);
        Assert.assertEquals(expected.size(),out.size());
        for (int i = 0; i < expected.size(); i++) {
            Assert.assertEquals(expected.get(i), out.get(i));
        }
    }

    @Test
    public void testConvertCodeBlock() {
        final List<String> inline = new ArrayList<String>();
        inline.add("Deep and important commentary");
        inline.add("<source>");
        inline.add("final double relativeAccuracy = 1.0e-12;");
        inline.add("final double absoluteAccuracy = 1.0e-8");
        inline.add("final int    maxOrder         = 5");
        inline.add("</source>");
        inline.add("And the saga continues");
        final List<String> expected = new ArrayList<String>();
        expected.add("Deep and important commentary");
        expected.add("");
        expected.add("    final double relativeAccuracy = 1.0e-12;");
        expected.add("    final double absoluteAccuracy = 1.0e-8");
        expected.add("    final int    maxOrder         = 5");
        expected.add("And the saga continues");
        final List<String> out = Utils.convertSourceBlocksAndTrim(inline);
        Assert.assertEquals(expected.size(),out.size());
        for (int i = 0; i < expected.size(); i++) {
            Assert.assertEquals(expected.get(i), out.get(i));
        }
    }

    @Test
    public void testConvertInlineCode() {
        final List<String> inline = new ArrayList<String>();
        inline.add("Deep and important commentary");
        inline.add("  Now this is deep:  <source>final double relativeAccuracy = 1.0e-12;</source>  ");
        inline.add("And the saga continues");
        final List<String> expected = new ArrayList<String>();
        expected.add("Deep and important commentary");
        expected.add("Now this is deep:  `final double relativeAccuracy = 1.0e-12;`");
        expected.add("And the saga continues");
        final List<String> out = Utils.convertSourceBlocksAndTrim(inline);
        Assert.assertEquals(expected.size(),out.size());
        for (int i = 0; i < expected.size(); i++) {
            Assert.assertEquals(expected.get(i), out.get(i));
        }
    }

    @Test
    public void testStripXML() {
        final List<String> inline = new ArrayList<String>();
        inline.add(" <?xml version=1.0?>");
        inline.add("<!--");
        inline.add("Licensed to the Apache Software Foundation (ASF) under one or more");
        inline.add("blah");
        inline.add("blah");
        inline.add("<document>");
        inline.add("<body>");
        inline.add("<section>");     // 0
        inline.add("<subsection>");  // 1
        inline.add("deep stuff");    // 2
        inline.add("more deepness"); // 3
        inline.add("</subsection>"); // 4
        inline.add("</section>");    // 5
        inline.add("</body>");
        inline.add("</document>");
        final List<String> out = Utils.stripXML(inline);
        Assert.assertEquals(6, out.size());
        for (int i = 0; i < out.size(); i++) {
            Assert.assertEquals(inline.get(i + 7), out.get(i));
        }
    }

    @Test
    public void testConvertSingleLineLinks() {
        final List<String> inline = new ArrayList<String>();
        inline.add("Deep, important stuff");
        inline.add("a href but no good");
        inline.add("<a href=" + Utils.QUOTE + "../good/stuff.html" + Utils.QUOTE + ">Foo</a>");
        inline.add("blah");
        final List<String> out = Utils.convertLinks(inline);
        Assert.assertEquals(inline.size(), out.size());
        for (int i = 0; i < out.size(); i++) {
            if (i == 2) {
                Assert.assertEquals("[Foo](../good/stuff.html)", out.get(i));
            } else {
                Assert.assertEquals(inline.get(i), out.get(i));
            }
        }
    }

    @Test
    public void testConvertMultiLineLinks() {
        final List<String> inline = new ArrayList<String>();
        inline.add("Deep, important stuff");
        inline.add("a href but no good");
        inline.add("<a href=" + Utils.QUOTE + "../good/stuff.html" + Utils.QUOTE + ">");
        inline.add("All about the Foo</a> plus some Bar.");
        final List<String> out = Utils.convertLinks(inline);
        Assert.assertEquals(inline.size(), out.size());
        for (int i = 0; i < out.size(); i++) {
            if (i < 2) {
                Assert.assertEquals(inline.get(i), out.get(i));
            }
        }
        Assert.assertEquals("[All about the Foo](../good/stuff.html)", out.get(2));
        Assert.assertEquals(" plus some Bar.", out.get(3));
    }

    @Test
    public void testConvertMultipleLinksPerLine() {
        final List<String> inline = new ArrayList<String>();
        inline.add("Deep, important stuff");
        inline.add("a href but no good");
        inline.add("<a href=" + Utils.QUOTE + "../good/stuff.html" + Utils.QUOTE + ">");
        inline.add("All about the Foo</a> plus some Bar. <a href=" + Utils.QUOTE +
                   "../bad/stuff.html" + Utils.QUOTE + ">Bad Stuff</a>");
        inline.add("Now the really bad part.");
        inline.add("<a href=" + Utils.QUOTE + "../good.html" + Utils.QUOTE + ">Deep Link</a>" +
                   "<a href=" + Utils.QUOTE + "../good2.html" + Utils.QUOTE + ">Deep Link2</a>");
        final List<String> out = Utils.convertLinks(inline);
        Assert.assertEquals(inline.size(), out.size());
        for (int i = 0; i < out.size(); i++) {
            if (i < 2) {
                Assert.assertEquals(inline.get(i), out.get(i));
            }
        }
        Assert.assertEquals("[All about the Foo](../good/stuff.html)", out.get(2));
        Assert.assertEquals(" plus some Bar. [Bad Stuff](../bad/stuff.html)", out.get(3));
        Assert.assertEquals("Now the really bad part.", out.get(4));
        Assert.assertEquals("[Deep Link](../good.html)[Deep Link2](../good2.html)", out.get(5));
    }

    @Test
    public void testConvertHeadings() {
        final List<String> inline = new ArrayList<String>();
        inline.add("<section name=" + Utils.QUOTE + "4 Numerical Analysis" + Utils.QUOTE + ">");
        inline.add("<subsection name=" + Utils.QUOTE + "4.1 Overview" + Utils.QUOTE +
                   " href=" + Utils.QUOTE + "overview" + Utils.QUOTE + ">");
        inline.add("<p>");
        inline.add("Analysis is fun.");
        inline.add("I mean, really fun.");
        inline.add("</subsection>");
        inline.add("<subsection name=" + Utils.QUOTE + "4.2 Meaty" + Utils.QUOTE +
                   " href=" + Utils.QUOTE + "overview" + Utils.QUOTE + ">");
        inline.add("Now, the meat.");
        inline.add("</subsection>");
        inline.add("</section>");
        final List<String> out = Utils.convertHeadings(inline);
        Assert.assertEquals(inline.size() - 1, out.size());
        Assert.assertEquals("# 4 Numerical Analysis", out.get(0));
        Assert.assertEquals("## 4.1 Overview", out.get(1));
        for (int i = 2; i < 5; i++) {
            Assert.assertEquals(inline.get(i), out.get(i));
        }
        Assert.assertEquals("", out.get(5));
        Assert.assertEquals("## 4.2 Meaty", out.get(6));
        Assert.assertEquals(inline.get(7), out.get(7));
        Assert.assertEquals("", out.get(8));
    }

    @Test
    public void testConvertCodeElements() {
        final List<String> inline = new ArrayList<String>();
        inline.add("Blah <code>x++</code>");
        inline.add("Love to code.");
        inline.add("<code>y++</code> to <code>m</code> and <code>n</code>");
        final List<String> out = Utils.convertCodeElements(inline);
        Assert.assertEquals("Blah `x++`", out.get(0));
        Assert.assertEquals(inline.get(1), out.get(1));
        Assert.assertEquals("`y++` to `m` and `n`", out.get(2));
    }

    @Test
    public void testConvertStandardParagraph() {
        final List<String> inline = new ArrayList<String>();
        inline.add("  <p>");
        inline.add("     Math is fun.");
        inline.add("  </p>");
        final List<String> out = Utils.convertParagraphs(inline);
        Assert.assertEquals(2, out.size());
        Assert.assertEquals(inline.get(1), out.get(0));
        Assert.assertEquals("", out.get(1));
    }

    @Test
    public void testConvertParagraphInlineStart() {
        final List<String> inline = new ArrayList<String>();
        inline.add("  <p>Math is fun.");
        inline.add("     Yes, really fun.");
        inline.add("  </p>");
        final List<String> out = Utils.convertParagraphs(inline);
        Assert.assertEquals(3, out.size());
        Assert.assertEquals("Math is fun.", out.get(0));
        Assert.assertEquals(inline.get(1), out.get(1));
        Assert.assertEquals("", out.get(2));
    }

    @Test
    public void testConvertParagraphInlineStartEnd() {
        final List<String> inline = new ArrayList<String>();
        inline.add("  <p>Math is fun.");
        inline.add("     Yes, really fun.</p>");
        final List<String> out = Utils.convertParagraphs(inline);
        Assert.assertEquals(3, out.size());
        Assert.assertEquals("Math is fun.", out.get(0));
        Assert.assertEquals("     Yes, really fun.", out.get(1));
        Assert.assertEquals("", out.get(2));
    }

    @Test
    public void testConvertParagraphOneLine() {
        final List<String> inline = new ArrayList<String>();
        inline.add("  <p>Math is fun.</p>");
        final List<String> out = Utils.convertParagraphs(inline);
        Assert.assertEquals(2, out.size());
        Assert.assertEquals("Math is fun.", out.get(0));
        Assert.assertEquals("", out.get(1));
    }

    @Test
    public void testConvertTableNoHeader() {
        final List<String> inline = new ArrayList<String>();
        inline.add("Deep stuff before table...");
        inline.add("  <table>");
        inline.add("    <tr><td>Hello</td><td>World</td></tr>");
        inline.add("    <tr><td>Life is</td><td> Good</td></tr>");
        inline.add("  </table>");
        inline.add("Deep stuff continues.");
        final List<String> out = Utils.convertTables(inline);
        Assert.assertEquals(5, out.size());
        Assert.assertEquals(inline.get(0), out.get(0));
        Assert.assertEquals(inline.get(5), out.get(4));
        Assert.assertEquals("| Hello | World |", out.get(1));
        Assert.assertEquals("| --- | --- |", out.get(2));
        Assert.assertEquals("| Life is | Good |", out.get(3));
    }

    @Test
    public void testConvertTableWithHeader() {
        final List<String> inline = new ArrayList<String>();
        inline.add("Deep stuff before table...");
        inline.add("  <table>");
        inline.add("    <tr><th>Hello</th><th>World</th></tr>");
        inline.add("    <tr><td>Life is</td><td> Good</td></tr>");
        inline.add("  </table>");
        inline.add("Deep stuff continues.");
        final List<String> out = Utils.convertTables(inline);
        Assert.assertEquals(5, out.size());
        Assert.assertEquals(inline.get(0), out.get(0));
        Assert.assertEquals(inline.get(5), out.get(4));
        Assert.assertEquals("| Hello | World |", out.get(1));
        Assert.assertEquals("| --- | --- |", out.get(2));
        Assert.assertEquals("| Life is | Good |", out.get(3));
    }

    @Test
    public void testConvertULSimple() {
        final List<String> inline = new ArrayList<String>();
        inline.add("Deep stuff before list...");
        inline.add("  <ul>");
        inline.add("    <li>one</li>");
        inline.add("    <li>two</li>");
        inline.add("  </ul>");
        inline.add("Deep stuff continues.");
        final List<String> out = Utils.convertLists(inline);
        Assert.assertEquals(5,out.size());
        Assert.assertEquals(inline.get(0), out.get(0));
        Assert.assertEquals(inline.get(inline.size() - 1), out.get(out.size() - 1));
        Assert.assertEquals("* one", out.get(1));
        Assert.assertEquals("* two", out.get(2));
        Assert.assertEquals(" ", out.get(3));
    }

    @Test
    public void testConvertULEmbeddedMarkup() {
        final List<String> inline = new ArrayList<String>();
        inline.add("Deep stuff before list...");
        inline.add("  <ul>");
        inline.add("    <li><p>one</p></li>");
        inline.add("    <li>two</li>");
        inline.add("  </ul>");
        inline.add("Deep stuff continues.");
        final List<String> out = Utils.convertLists(inline);
        Assert.assertEquals(5,out.size());
        Assert.assertEquals(inline.get(0), out.get(0));
        Assert.assertEquals(inline.get(inline.size() - 1), out.get(out.size() - 1));
        Assert.assertEquals("* <p>one</p>", out.get(1));
        Assert.assertEquals("* two", out.get(2));
        Assert.assertEquals(" ", out.get(3));
    }

    @Test
    public void testConvertOLSimple() {
        final List<String> inline = new ArrayList<String>();
        inline.add("Deep stuff before list...");
        inline.add("  <ol>");
        inline.add("    <li>one</li>");
        inline.add("    <li>two</li>");
        inline.add("  </ol>");
        inline.add("Deep stuff continues.");
        final List<String> out = Utils.convertLists(inline);
        Assert.assertEquals(5,out.size());
        Assert.assertEquals(inline.get(0), out.get(0));
        Assert.assertEquals(inline.get(inline.size() - 1), out.get(out.size() - 1));
        Assert.assertEquals("1. one", out.get(1));
        Assert.assertEquals("1. two", out.get(2));
        Assert.assertEquals(" ", out.get(3));
    }

    @Test
    public void testConvertNestedLists() {
        final List<String> inline = new ArrayList<String>();
        inline.add("Deep stuff before list...");
        inline.add("  <ol>");
        inline.add("    <li>zero</li>");
        inline.add("    <li><ul>");
        inline.add("            <li>one</li>");
        inline.add("            <li>two</li>");
        inline.add("       </ul>");
        inline.add("    <li>three</li>");
        inline.add("    <li>four</li>");
        inline.add("  </ol>");
        inline.add("Deep stuff continues.");
        final List<String> out = Utils.convertLists(inline);
        final List<String> expected = new ArrayList<String>();
        expected.add("Deep stuff before list...");
        expected.add("1. zero");
        expected.add("    * one");
        expected.add("    * two");
        expected.add("1. three");
        expected.add("1. four");
        expected.add(" ");
        expected.add("Deep stuff continues.");
        Assert.assertEquals(expected, out);
    }
}