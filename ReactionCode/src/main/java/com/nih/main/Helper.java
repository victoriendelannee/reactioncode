/*
 * Copyright (C) 2007-2018 Syed Asad Rahman <asad @ ebi.ac.uk>.
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA 02110-1301  USA
 */
package com.nih.main;

import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintWriter;
import static java.lang.System.getProperty;
import static java.lang.System.out;

import java.util.Map;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;

import com.nih.tools.ColouredSystemOutPrintln;

/**
 * @contact Syed Asad Rahman, EMBL-EBI, Cambridge, UK.
 * @author Syed Asad Rahman <asad @ ebi.ac.uk>
 */
class Helper {

    static final String NEW_LINE = getProperty("line.separator");

    protected static void getHeader() {
        StringBuilder sb = new StringBuilder();

        sb.append("**************************************************");
        sb.append(NEW_LINE);
        sb.append("Pseudo-Molecule & ReactionCode");
        sb.append(NEW_LINE);
        sb.append(NEW_LINE);
        sb.append("Author: Victorien Delannee");
        sb.append(NEW_LINE);
        sb.append("e-mail: victorien.delannee@nih.gov");
        sb.append(NEW_LINE);
        sb.append("NCI-National Cancer Institute");
        sb.append(NEW_LINE);
        sb.append("Frederick, Maryland");
        sb.append(NEW_LINE);
        sb.append("USA");
        sb.append(NEW_LINE);
        sb.append(NEW_LINE);
        sb.append("Note: The copyright of this software belongs to the author");
        sb.append(NEW_LINE);
        sb.append("and NIH-National Institutes of Health");
        sb.append(NEW_LINE);
        sb.append(NEW_LINE);

        sb.append("Reference");
        sb.append(NEW_LINE);
        sb.append("Victorien Delannée and Marc Nicklaus");
        sb.append(NEW_LINE);
        sb.append("ReactionCode: a new versatile format for searching, analysis, classification, transform, and encoding/decoding of reactions");
        sb.append(NEW_LINE);
        sb.append("doi: 10.26434/chemrxiv.12058971");
        sb.append(NEW_LINE);
        sb.append("**************************************************");
        sb.append(NEW_LINE);
        
        out.println(sb.toString());
    }

    /**
     * WreactionWithLayoutite the preactionWithLayoutovided
     * numbereactionWithLayout of blank lineheaderString to the
     * preactionWithLayoutovided OutputStreactionWithLayouteam.
     *
     * @param numberBlankLines NumbereactionWithLayout of blank lineheaderString
     * to wreactionWithLayoutite.
     * @param out OutputStreactionWithLayouteam to which to
     * wreactionWithLayoutite the blank lineheaderString.
     */
    protected static void displayBlankLines(final int numberBlankLines, final OutputStream out) {
        try {
            for (int i = 0; i < numberBlankLines; ++i) {
                out.write(NEW_LINE.getBytes());
            }
        } catch (IOException ioEx) {
            for (int i = 0; i < numberBlankLines; ++i) {
                System.out.println();
            }
        }
    }

    /*
     System.out.println("-- USAGE --");
     printUsage(applicationName + " (Posix)", constructPosixOptions(), System.out);
     displayBlankLines(1, System.out);
     printUsage(applicationName + " (Gnu)", constructGnuOptions(), System.out);
     displayBlankLines(4, System.out);
     System.out.println("-- HELP --");
     printHelp(
     constructPosixOptions(), 80, "POSIX HELP", "End of POSIX Help",
     3, 5, true, System.out);
     displayBlankLines(1, System.out);
     printHelp(
     constructGnuOptions(), 80, "GNU HELP", "End of GNU Help",
     5, 3, true, System.out);
     */
    protected static void printHelp(final OutputStream out, final Options options) {
        final String commandLineSyntax = "java -jar ReactionCode-1.0.0.jar";
        try (PrintWriter writer = new PrintWriter(out)) {
            final HelpFormatter formatter = new HelpFormatter();
            displayBlankLines(2, out);
            formatter.printHelp(writer, 80, commandLineSyntax, "HELP",
                    options, 5, 3, "End of Helper Help", true);
            writer.flush();
            writer.close();
        }
    }

    protected static void printHelp(final Map<String, Options> optionsMap, final int printedRowWidth,
            final String header, final String footer, final int spacesBeforeOption,
            final int spacesBeforeOptionDescription, final boolean displayUsage, final OutputStream out) throws IOException {
        final String commandLineSyntax = "java -jar ReactionCode-1.0.0.jar";
        try (PrintWriter writer = new PrintWriter(out)) {
            final HelpFormatter helpFormatter = new HelpFormatter();
            //do not sort, comment next line to sort by alphabetic order
            helpFormatter.setOptionComparator(null);
            optionsMap.keySet().stream().map((headerString) -> {
                helpFormatter.printHelp(
                        writer,
                        printedRowWidth,
                        commandLineSyntax,
                        headerString,
                        optionsMap.get(headerString),
                        spacesBeforeOption,
                        spacesBeforeOptionDescription,
                        "End of Helper " + headerString + " Help",
                        displayUsage);
                return headerString;
            }).map((_item) -> {
                displayBlankLines(2, out);
                return _item;
            }).forEach((_item) -> {
                writer.flush();
            });
            writer.close();
        }
    }

    protected static void printUsageExamples() {
        StringBuilder sb = new StringBuilder();
        sb.append(NEW_LINE);
        sb.append("Encoder examples: ");
        sb.append(NEW_LINE);
        sb.append("java -jar ReactionCode-1.2.0.jar -q reaction.rxn -o /path/to/output/");
        sb.append(NEW_LINE);
        sb.append("java -jar ReactionCode-1.2.0.jar -q reaction.rdf");
        sb.append(NEW_LINE);
        sb.append("java -jar ReactionCode-1.2.0.jar -q \"[cH2:1]1([Br:21])[cH2:2][cH2:3][cH:4]([cH2:5][cH2:6]1)[cH:7]2[cH2:8][cH2:9][cH2:10]([Br:22])[cH2:11][cH2:12]2.[NH2:13][cH:14]3[cH2:15][cH2:16][cH2:17][cH:18]([cH2:19]3)[CH3:20]>>[cH2:1]1([NH2:13][cH:14]3[cH2:15][cH2:16][cH2:17][cH:18]([cH2:19]3)[CH3:20])[cH2:2][cH2:3][cH:4]([cH2:5][cH2:6]1)[cH:7]2[cH2:8][cH2:9][cH2:10]([NH2:13][cH:14]3[cH2:15][cH2:16][cH2:17][cH:18]([cH2:19]3)[CH3:20])[cH2:11][cH2:12]2\"");
        sb.append(NEW_LINE);
        sb.append(NEW_LINE);
        sb.append("Decoder examples: ");
        sb.append(NEW_LINE);
        sb.append("java -jar ReactionCode-1.2.0.jar -q reactionCodes.txt -r -s");
        sb.append(NEW_LINE);
        sb.append("java -jar ReactionCode-1.2.0.jar -q \"0:90E()[1]908(01GG)[1]708(10GG)[1]506(21GH)[1]506(12GJ)[1]|1:006(11GJ)[1]006(11GG)[1]006(11GG)[1]006(11GG)[1]/i06II|2:008(11GM)[1]008(11GM)[1]006(22GL)[1]006(11GL)[1]006(11GM)[1]/c00HH/i00JJ02HH|3:006(11GS)[1]006(11GU11GR)[1]/s00210212|4:008(11GU)[1]006(11GV)[1]/s01640364|5:006(11GW)[1]006(22GX)[1]006(11GX)[1]|A:010(11GI)[1]|B:008(2200)[1]008(2200)[1]006(1100)[1]|C:009(1103)[1]009(1103)[1]009(1103)[1]|\" -m -i -p \"decode\"");
        sb.append(NEW_LINE);
        sb.append(NEW_LINE);
        sb.append("Transformer example: ");
        sb.append(NEW_LINE);
        sb.append("java -jar ReactionCode-1.2.0.jar -q \"0:907()[1]906(01GG)[1]711(10GH)[1]|1:007(99GH)[1]006(99GH)[1]006(11GG)[1]|\" -t \"[NH2:1][CH:6]1[CH2:11][CH2:17][CH2:22][CH2:18][CH2:12]1.[Cl:2][c:3]1[n:4]([n:7][c:8]2[c:5]1[cH:10][cH:16][c:19]([cH:13]2)[F:23])-[c:9]3[cH:14][cH:20][c:24]([cH:21][cH:15]3)[Cl:25]\" -m ");
        sb.append(NEW_LINE);
        out.println(sb.toString());
    }

}
