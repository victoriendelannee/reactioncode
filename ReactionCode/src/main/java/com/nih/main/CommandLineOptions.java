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

import org.apache.commons.cli.Options;

public class CommandLineOptions {

    /**
     *
     */
    public CommandLineOptions() {
    }

    /**
     *
     * @return
     */
    protected Options createEncoderCodeOptions() {
        Options options = new Options();
        options.addOption("q", "query", true, "RXN/RDF/SMIRKS (file or String)");
        options.addOption("o", "output", true, "Output directory");
        options.addOption("p", "prefix", true, "Job prefix");
        options.addOption("f", "format", true, "Output format (JSON/CSV)");
        options.addOption("g", "getPseudoSmiles", false, "Create SMILES of the Pseudo-Molecules");
        options.addOption("i", "image", false, "Create image of the Pseudo-Molecules");
        options.addOption("s", "sdf", false, "Write SDF of the Pseudo-Molecules");
        options.addOption("n", "onlyPseudoMolecule", false, "Generate only the Pseudo-Molecules");
        options.addOption("a", "alternativeVersion", false, "Alternative version of the ReactionCode");
        options.addOption("r", "norepetiton", false, "To not set the atom repetition in the ReactionCode");
        options.addOption("e", "examples", false, "Show examples");
        options.addOption("h", "help", false, "Help page for command usage");

        return options;
    }
    
    /**
    *
    * @return
    */
   protected Options createDecoderCodeOptions() {
       Options options = new Options();
       options.addOption("q", "query", true, "Multi-lines ReactionCode file (one ReactionCode per line) or ReactionCode String");
       options.addOption("o", "output", true, "Output directory");
       options.addOption("p", "prefix", true, "Job prefix");
       options.addOption("i", "image", false, "Create image of the Reaction");
       options.addOption("r", "RDF", false, "Write RDF of the Reaction");
       options.addOption("x", "RXN", false, "Write RXN of the Reaction");
       options.addOption("m", "SMIRKS", false, "Write SMIRKS (reaction SMARTS)");
       options.addOption("s", "SMILES", false, "Write Reaction SMILES");
       options.addOption("e", "examples", false, "Show examples");
       options.addOption("h", "help", false, "Help page for command usage");

       return options;
   }
   
   /**
   *
   * @return
   */
  protected Options createTransformerCodeOptions() {
      Options options = new Options();
      options.addOption("q", "query", true, "Multi-lines ReactionCode file (one ReactionCode per line) or ReactionCode String");
      options.addOption("t", "transformer", true, "Reactants (SMILES)");
      options.addOption("o", "output", true, "Output directory");
      options.addOption("p", "prefix", true, "Job prefix");
      options.addOption("a", "all", false, "generate all possible reactions");
      options.addOption("u", "unique", false, "generate one possible reaction (best match)");
      options.addOption("z", "stoichiometry", true, "stoichiometry (number of times a same reactant can react on different site)");
      options.addOption("c", "checkValence", false, "Check the valence of the atom in reaction center to filter candidates");
      options.addOption("i", "image", false, "Create image of the Reaction");
      options.addOption("r", "RDF", false, "Write RDF of the Reaction");
      options.addOption("x", "RXN", false, "Write RXN of the Reaction");
      options.addOption("m", "SMIRKS", false, "Write SMIRKS (reaction SMARTS)");
      options.addOption("s", "SMILES", false, "Write Reaction SMILES");
      options.addOption("e", "examples", false, "Show examples");
      options.addOption("h", "help", false, "Help page for command usage");

      return options;
  }
}
