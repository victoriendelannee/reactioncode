/*
 * Copyright (C) 2007-2018 Syed Asad Rahman <asad @ ebi.ac.uk>.
 * (modified by victorien delannee)
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
package com.nih.parser;

import static java.lang.String.format;
import static java.lang.System.exit;
import static java.util.logging.Level.INFO;
import static java.util.logging.Level.SEVERE;
import static java.util.logging.Level.WARNING;
import static org.openscience.cdk.io.IChemObjectReader.Mode.RELAXED;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.regex.Pattern;

import org.openscience.cdk.AtomContainer;
import org.openscience.cdk.Reaction;
import org.openscience.cdk.ReactionSet;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IReaction;
import org.openscience.cdk.interfaces.IReactionSet;
import org.openscience.cdk.io.CMLReader;
import org.openscience.cdk.io.Mol2Reader;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.smiles2.SmilesParser2;
import org.openscience.cdk.tools.ILoggingTool;
import org.openscience.cdk.tools.LoggingToolFactory;

import com.nih.codes.DecodeReactionCode;
import com.nih.tools.ColouredSystemOutPrintln;
import com.nih.writer.MDLV2000RDFWriter;
import com.opencsv.CSVParser;
import com.opencsv.CSVParserBuilder;
import com.opencsv.CSVReader;
import com.opencsv.CSVReaderBuilder;


public class ChemicalFormatParser {

	private static final ILoggingTool LOGGER
	= LoggingToolFactory.createLoggingTool(ChemicalFormatParser.class);

	List<File> files = new ArrayList<File>();
	int reactionCount;
	private boolean kekulize;
	private boolean radicalize = false;

	/**
	 * Detect the format of the input. It can be a string like SMILES, SMIRKS or ReactionCode and also a file
	 * like RDF, RXN, a file containing one ReactionCode per line, or one SMIRKS per line, or a CSV with
	 * a column containg the ReactionCOdes
	 * @param fileName
	 * @return
	 * @throws IOException
	 */
	public String formatDetector(String fileName) throws IOException {
		String format = null;
		File file = new File(fileName); 
		if (file.isFile()) {
			BufferedReader br = new BufferedReader(new FileReader(file)); 
			String str; 
			while ((str = br.readLine()) != null) {
				if (!str.isEmpty()) {
					if (str.contains("$RXN")) {
						format = "RXN_FILE";
						break;
					}
					else if (str.contains("$RDFILE")) {
						format = "RDF_FILE";
						break;
					}
					else if (str.contains(">>")) {
						format = "SMIRKS_FILE";
						break;
					}
					else if (str.contains("|") && str.contains(":") && str.contains("(")) {
						format = "REACTIONCODE_FILE";
					}
					else if (Pattern.compile(" ").matcher(str).find()) {
						String[] split;
						String separator = "SPACE";
						split = str.split(" ");
						for (int i = 0; i < split.length; i++) {
							String column = split[i];
							if (column.contains("|") && column.contains(":") && column.contains("(")) {
								format = "REACTIONCODE_CSV_WITHOUTHEADER_"+i+"_"+separator;
								break;
							}
						}
						break;
					}
					else if (Pattern.compile("\t").matcher(str).find()) {
						String[] split;
						String separator = "TAB";
						split = str.split("\t");
						for (int i = 0; i < split.length; i++) {
							String column = split[i];
							if (column.contains("|") && column.contains(":") && column.contains("(")) {
								format = "REACTIONCODE_CSV_WITHOUTHEADER_"+i+"_"+separator;
								break;
							}
						}
						break;
					}
					else if (Pattern.compile(",").matcher(str).find()) {
						String[] split;
						String separator = ",";
						split = str.split(",");
						for (int i = 0; i < split.length; i++) {
							String column = split[i];
							if (column.contains("|") && column.contains(":") && column.contains("(")) {
								format = "REACTIONCODE_CSV_WITHOUTHEADER_"+i+"_"+separator;
								break;
							}
						}
						break;
					}
					else {
						String nextline = br.readLine();
						String[] split;
						String separator = "";
						if (Pattern.compile(" ").matcher(nextline).find())  {
							split = nextline.split(" ");
							separator = "SPACE";
						}
						else if (Pattern.compile("\t").matcher(nextline).find()) {
							split = nextline.split("\t");
							separator = "TAB";
						}
						else if (nextline.contains(",")) {
							split = nextline.split(",");
							separator = ",";
						}
						else {
							break;
						}
						for (int i = 0; i < split.length; i++) {
							String column = split[i];
							if (column.contains("|") && column.contains(":") && column.contains("(")) {
								format = "REACTIONCODE_CSV_WITHHEADER_"+i+"_"+separator;
								break;
							}
						}
						break;
					}
				}
			}
			br.close();
		}
		else {
			if (fileName.contains(">>")) {
				format = "SMIRKS_TEXT";
			}
			else if (fileName.contains("|") && fileName.contains(":") && fileName.contains("(")) {
				format = "REACTIONCODE_TEXT";
			}
		}
		return format; 
	}

	/**
	 * RDF parser
	 * @param fileName
	 * @return
	 * @throws CDKException
	 * @throws IOException
	 */
	public IReactionSet parseRDF(String fileName) throws CDKException, IOException {
		InputStream ins = new FileInputStream(fileName);
		//System.out.println(IOUtils.toString(ins));
		MDLReactionsFileReader reader = new MDLReactionsFileReader(ins);
		reader.setAromatize(true);
		IReactionSet reactions = (IReactionSet) reader.read(new ReactionSet());
		reader.close();
		return reactions;
	}
	
	/**
	 * RDF parser which split by chunk of 5000
	 * @param fileName
	 * @param split
	 * @return
	 * @throws CDKException
	 * @throws IOException
	 */
	public IReactionSet parseRDF(String fileName, boolean split) throws CDKException, IOException {
		InputStream ins = new FileInputStream(fileName);
		//System.out.println(IOUtils.toString(ins));
		MDLReactionsFileReader reader = new MDLReactionsFileReader(ins);
		reader.setAromatize(true);
		if (split == true) {
			reader.setSplit(true);
			reader.setFilename(fileName);
		}
		IReactionSet reactions = (IReactionSet) reader.read(new ReactionSet());
		if (split == true)
			files = reader.getFiles();
		reactionCount = reader.getReactionCount();
		reader.close();
		return reactions;
	}

	/**
	 * RXN parser
	 * @param fileNames
	 * @return
	 * @throws CDKException
	 * @throws IOException
	 */
	public IReactionSet parseRXN(String fileNames) throws CDKException, IOException {
		IReactionSet reactionsSet = DefaultChemObjectBuilder.getInstance().newInstance(IReactionSet.class);
		String[] files = fileNames.split(";");

		for (String file : files) {
			InputStream ins = new FileInputStream(file);
			//System.out.println(IOUtils.toString(ins));
			MDLReactionsFileReader reader = new MDLReactionsFileReader(ins);
			reader.setAromatize(true);
			IReaction reaction = (IReaction) reader.read(new Reaction());
			reader.close();

			reactionsSet.addReaction(reaction);
		}
		reactionCount = reactionsSet.getReactionCount();
		return reactionsSet;
	}

	/**
	 * @param input
	 * @return
	 * @throws CDKException
	 * @throws IOException
	 */
	public IReaction parseCML(String input) throws CDKException, IOException {
		File f = new File(input);
		if (!f.isFile()) {
			LOGGER.warn(WARNING, format("CML file not found! " + f.getName()));
			exit(1);
		}
		String[] split = f.getName().split(".cml");
		CMLReader cmlReader = new CMLReader(new FileInputStream(input));
		AtomContainer ac = cmlReader.read(new AtomContainer());
		cmlReader.close();
		IReaction r = new Reaction();
		r.addReactant(ac, 1.0);
		r.addProduct(ac, 1.0);
		r.setID(split[0]);
		return r;
	}

	/**
	 * @param r
	 * @return
	 * @throws CDKException
	 */
	protected IReaction convertRoundTripRXNSMILES(IReaction r) throws CDKException {
		final SmilesGenerator sg = new SmilesGenerator(
				SmiFlavor.AtomAtomMap
				| SmiFlavor.UseAromaticSymbols
				| SmiFlavor.Stereo);
		String createSmilesFromReaction = sg.create(r);
		final SmilesParser2 smilesParser = new SmilesParser2(DefaultChemObjectBuilder.getInstance());
		IReaction parseReactionSmiles = smilesParser.parseReactionSmiles(createSmilesFromReaction);
		parseReactionSmiles.setID(r.getID());
		for (int i = 0; i < r.getReactantCount(); i++) {
			parseReactionSmiles.getReactants().getAtomContainer(i).setID(r.getReactants().getAtomContainer(i).getID());
		}
		for (int i = 0; i < r.getProductCount(); i++) {
			parseReactionSmiles.getProducts().getAtomContainer(i).setID(r.getProducts().getAtomContainer(i).getID());
		}
		reactionCount = 1;
		return parseReactionSmiles;
	}

	/**
	 * Reaction SMILES and SMIRKS parser
	 * @param file
	 * @param id
	 * @return
	 * @throws IOException
	 */
	public IReactionSet parseReactionSMILES(String file, int id) throws IOException {
		//List<String> reactionSmiles = new ArrayList<String>();
		BufferedReader br = new BufferedReader(new FileReader(file)); 

		String smiles; 
		SmilesParser2 sp = new SmilesParser2(DefaultChemObjectBuilder.getInstance());
		sp.kekulise(kekulize);
		sp.radicalize(radicalize);
		IReactionSet reactions = DefaultChemObjectBuilder.getInstance().newInstance(IReactionSet.class); 
		int smilesIndex = id + 1;
		while ((smiles = br.readLine()) != null) {
			if (smiles.contains(">>")) {
				smiles = smiles.replaceAll("[^\\x00-\\x7F]", "");
				try {
					IReaction parseReactionSmile = sp.parseReactionSmiles(smiles);
					try {
//						LOGGER.error(INFO, "Annotating Reaction " + "smiles");
//						if (reactionSmiles.size() > 1) {
							parseReactionSmile.setID("smiles_" + smilesIndex);
//						} else {
//							parseReactionSmile.setID("smiles");
//						}
						Map<Object,Object> rdFields = new HashMap<Object,Object>();
						rdFields.put("SMIRKS", smiles);
						parseReactionSmile.addProperties(rdFields);
						reactions.addReaction(parseReactionSmile);
					} catch (Exception ex) {
						LOGGER.error(SEVERE, null, ex);
					}
				} catch (InvalidSmilesException ex) {
					LOGGER.error(SEVERE, null, ex);
				}
				smilesIndex++;
			}
		} 
		br.close();

		
		reactionCount = reactions.getReactionCount();
		return reactions;
	}
	
	/**
	 * Reaction SMILES and SMIRKS parser
	 * @param file
	 * @param id
	 * @param split
	 * @return
	 * @throws IOException
	 */
	public IReactionSet parseReactionSMILES(String file, int id, boolean split) throws IOException {
		Path path = Paths.get(file);
		long lineCount = Files.lines(path).count();
		
		File directory = new File(new File(file).getParent(), "tempRDFiles/");
		if (lineCount > 10000) {
			if (!directory.exists()){
		        directory.mkdir();
		    }
		}
		
		FileWriter writer = null;
		
		//List<String> reactionSmiles = new ArrayList<String>();
		BufferedReader br = new BufferedReader(new FileReader(file)); 

		String smiles; 
		SmilesParser2 sp = new SmilesParser2(DefaultChemObjectBuilder.getInstance());
		sp.kekulise(kekulize);
		sp.radicalize(radicalize);
		IReactionSet reactions = DefaultChemObjectBuilder.getInstance().newInstance(IReactionSet.class); 
		int smilesIndex = id + 1;
		int tempReactionCounter = 0;
		while ((smiles = br.readLine()) != null) {
			if (smiles.contains(">>")) {
				smiles = smiles.replaceAll("[^\\x00-\\x7F]", "");
				if (split == true && reactionCount >= 10000) {
					if (tempReactionCounter == 10000) {
						if (writer != null)
							writer.close();
						File temp = new File(directory.toString() + "/chunk_" + files.size() + ".txt");
						files.add(temp);
						writer = new FileWriter(temp);
						writer.write(smiles + System.lineSeparator());
						tempReactionCounter = 0;
						System.out.println(ColouredSystemOutPrintln.ANSI_GREEN + reactionCount + 
								" reactions are counted and a temp file has been created" + 
								ColouredSystemOutPrintln.ANSI_RESET);
					}
					else {
						writer.write(smiles + System.lineSeparator());
					}
				}
				else if (reactionCount < 10000 || split == false) {
					try {
						IReaction parseReactionSmile = sp.parseReactionSmiles(smiles);
						try {
							//LOGGER.error(INFO, "Annotating Reaction " + "smiles");
							//if (reactionSmiles.size() > 1) {
							parseReactionSmile.setID("smiles_" + smilesIndex);
							//} else {
							//parseReactionSmile.setID("smiles");
							//}
							Map<Object,Object> rdFields = new HashMap<Object,Object>();
							rdFields.put("SMIRKS", smiles);
							parseReactionSmile.addProperties(rdFields);
							reactions.addReaction(parseReactionSmile);
						} catch (Exception ex) {
							LOGGER.error(SEVERE, null, ex);
						}
					} catch (InvalidSmilesException ex) {
						LOGGER.error(SEVERE, null, ex);
					}
					smilesIndex++;
				}
				reactionCount++;
				tempReactionCounter++;
			}
		} 
		br.close();
		
		if (writer != null)
			writer.close();

		return reactions;
	}

	
	/**
	 * Reaction SMILES and SMIRKS parser
	 * @param file
	 * @param id
	 * @param split
	 * @param kekule
	 * @return
	 * @throws IOException
	 */
//	public IReactionSet parseReactionSMILES(String file, int id, boolean split, boolean kekule) throws IOException {
//		Path path = Paths.get(file);
//		long lineCount = Files.lines(path).count();
//		
//		//List<String> reactionSmiles = new ArrayList<String>();
//		BufferedReader br = new BufferedReader(new FileReader(file)); 
//
//		String smiles; 
//		SmilesParser sp = new SmilesParser(getInstance());
//		sp.kekulise(kekule);
//		IReactionSet reactions = DefaultChemObjectBuilder.getInstance().newInstance(IReactionSet.class); 
//		int smilesIndex = id + 1;
//		while ((smiles = br.readLine()) != null) {
//			if (smiles.contains(">>")) {
//				smiles = smiles.replaceAll("[^\\x00-\\x7F]", "");
//				try {
//					if (split == true && reactions.getReactionCount() == 10000) {
//						writeRDFChunk(reactions, file);
//						reactions.removeAllReactions();
//						System.out.println(ColouredSystemOutPrintln.ANSI_GREEN + reactionCount + 
//								" reactions are parsed and a tempory RD file has been created" + 
//								ColouredSystemOutPrintln.ANSI_RESET);
//					}
//					IReaction parseReactionSmile = sp.parseReactionSmiles(smiles);
//					reactionCount++;
//					try {
////						LOGGER.error(INFO, "Annotating Reaction " + "smiles");
////						if (reactionSmiles.size() > 1) {
//							parseReactionSmile.setID("smiles_" + smilesIndex);
////						} else {
////							parseReactionSmile.setID("smiles");
////						}
//						Map<Object,Object> rdFields = new HashMap<Object,Object>();
//						rdFields.put("SMIRKS", smiles);
//						parseReactionSmile.addProperties(rdFields);
//						reactions.addReaction(parseReactionSmile);
//					} catch (Exception ex) {
//						LOGGER.error(SEVERE, null, ex);
//					}
//				} catch (InvalidSmilesException ex) {
//					LOGGER.error(SEVERE, null, ex);
//				}
//				smilesIndex++;
//			}
//		} 
//		br.close();
//
//		return reactions;
//	}
	

	/**
	 * Reaction SMILES and SMIRKS parser
	 * @param reactionSmiles
	 * @return
	 */
	public IReactionSet parseReactionSMILES(String reactionSmiles) {
		SmilesParser2 sp = new SmilesParser2(DefaultChemObjectBuilder.getInstance());
		sp.kekulise(kekulize);
		sp.radicalize(radicalize);
		String[] smiles = reactionSmiles.split("\\s+");
		IReactionSet reactions = DefaultChemObjectBuilder.getInstance().newInstance(IReactionSet.class); 
		int smilesIndex = 1;
		for (String s : smiles) {
			try {
				IReaction parseReactionSmile = sp.parseReactionSmiles(s);
				try {
//					LOGGER.error(INFO, "Annotating Reaction " + "smiles");
					if (smiles.length > 1) {
						parseReactionSmile.setID("smiles_" + smilesIndex);
					} else {
						parseReactionSmile.setID("smiles");
					}
					reactions.addReaction(parseReactionSmile);
				} catch (Exception ex) {
					LOGGER.error(SEVERE, null, ex);
				}
			} catch (InvalidSmilesException ex) {
				LOGGER.error(SEVERE, null, ex);
			}
			smilesIndex++;
		}
		reactionCount = reactions.getReactionCount();
		return reactions;
	}

	/**
	 * Parse a smiles containing multiple species (for example the left part of a reaction SMILES)
	 * F.N#Cc(c)c
	 * Do not kekule, if the SMILES/SMARTS is a part of the molecule (ex non complete aromatic ring)
	 * @param siles
	 * @return
	 * @throws InvalidSmilesException 
	 */
	public IAtomContainerSet parseSMILES(String smiles, boolean kekule) throws InvalidSmilesException {
		SmilesParser2 sp = new SmilesParser2(DefaultChemObjectBuilder.getInstance());
    	sp.kekulise(kekule);

    	String[] smis = smiles.split("\\.");
		IAtomContainerSet set = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);
    	for (String reactant : smis) {
    		IAtomContainer rea = sp.parseSmiles(reactant);
    		set.addAtomContainer(rea);
    	}
		return set;
	}
	
	/**
	 * Reaction SMILES and SMARTS parser
	 * @param smiles
	 * @return
	 */
	public IReaction parseSMILES(String smiles) {
		SmilesParser2 sp = new SmilesParser2(DefaultChemObjectBuilder.getInstance());
		try {
			IAtomContainer mol = sp.parseSmiles(smiles);
			try {
				IReaction parseReactionSmiles = DefaultChemObjectBuilder.getInstance().newInstance(IReaction.class);
				parseReactionSmiles.addReactant(mol, 1.0);
//				LOGGER.error(INFO, "Annotating Reaction " + "smiles");
				parseReactionSmiles.setID("smiles");
				reactionCount = 1;
				return parseReactionSmiles;
			} catch (IllegalArgumentException ex) {
				LOGGER.error(SEVERE, null, ex);
			}
		} catch (InvalidSmilesException ex) {
			LOGGER.error(SEVERE, null, ex);
		}
		return null;
	}

	/**
	 * MOL2 file parser
	 * @param input
	 * @return
	 * @throws CDKException
	 * @throws IOException
	 */
	public IReaction parseMOL2(String input) throws CDKException, IOException {
		File f = new File(input);
		if (!f.isFile()) {
			LOGGER.error(WARNING, format("Mol2 file not found! " + f.getName()));
			exit(1);
		}

		String[] split = f.getName().split(".mol");
		MDLV2000Reader mdlV2000Reader = new MDLV2000Reader(
				new FileReader(input), RELAXED);
		AtomContainer ac = mdlV2000Reader.read(new AtomContainer());
		mdlV2000Reader.close();
		IReaction r = new Reaction();
		r.addReactant(ac, 1.0);
		r.addProduct(ac, 1.0);
		r.setID(split[0]);
		return r;
	}

	/**
	 *SDF parser
	 * @param input
	 * @return
	 * @throws CDKException
	 * @throws IOException
	 */
	public IReaction parseSDF(String input) throws CDKException, IOException {
		File f = new File(input);
		if (!f.isFile()) {
			LOGGER.error(WARNING, format("SDF file not found! " + f.getName()));
			exit(1);
		}
		String[] split = f.getName().split(".sdf");
		Mol2Reader mol2Reader = new Mol2Reader(new FileReader(input));
		AtomContainer ac = mol2Reader.read(new AtomContainer());
		mol2Reader.close();
		IReaction r = new Reaction();
		r.addReactant(ac, 1.0);
		r.addProduct(ac, 1.0);
		r.setID(split[0]);
		return r;
	}

	/**
	 * ReactionCode parser (String)
	 * @param fileName
	 * @return
	 * @throws CDKException
	 * @throws IOException
	 * @throws CloneNotSupportedException
	 */
	public IReactionSet parseReactionCode(String fileName) throws CDKException, IOException, CloneNotSupportedException {
		List<String> reactionCodes = new ArrayList<String>();
		if (new File(fileName).isFile()) {
			BufferedReader br = new BufferedReader(new FileReader(fileName)); 
			String line; 
			String reactionCode = ""; 
			while ((line = br.readLine()) != null) {
				if (line.contains("\\:") && line.charAt(0) == '0') {
					if (reactionCode.length() > 0) 
						reactionCodes.add(reactionCode); 
					reactionCode = line;
				}
				else if (line.contains("\\:") && line.charAt(0) != '0') {
					reactionCode += line;
				}
			} 
			br.close();
		}
		DecodeReactionCode decoder = new DecodeReactionCode();
		IReactionSet reactionsSet = decoder.makeReactions(reactionCodes);
		
		return reactionsSet;
	}

	/**
	 * ReactionCode parser. Format has to be specified (REACTIONCODE_FILE, REACTIONCODE_TEXT, REACTIONCODE_CSV_)
	 * USE FORMAT DETECTOR IF NOT SURE
	 * @param fileName
	 * @param format
	 * @return
	 * @throws CDKException
	 * @throws IOException
	 * @throws CloneNotSupportedException
	 */
	public IReactionSet parseReactionCode(String fileName, String format) throws CDKException, IOException, CloneNotSupportedException {
		List<String> reactionCodes = new ArrayList<String>();
		
		if (format.equals("REACTIONCODE_FILE")) {
			BufferedReader br = new BufferedReader(new FileReader(fileName)); 
			String line; 
			String reactionCode = ""; 
			while ((line = br.readLine()) != null) {
				if (line.contains(":") && line.charAt(0) == '0') {
					if (reactionCode.length() > 0) 
						reactionCodes.add(reactionCode); 
					reactionCode = line;
				}
				else if (line.contains(":") && line.charAt(0) != '0') {
					reactionCode += line;
				}
			}
			if (reactionCode.length() > 0) 
				reactionCodes.add(reactionCode); 
			br.close();
		}
		else if (format.equals("REACTIONCODE_TEXT")) {
			reactionCodes.add(fileName);
		}
		else if (format.contains("REACTIONCODE_CSV_")) {
			CSVReader reader = null;
	        try {
	        	int column = Integer.parseInt(format.split("_")[3]);
	        	char separator = ','; 
	        	
	        	if (format.contains("SPACE")) {
	        		separator = ' ';
	        	}
	        	else if (format.contains("TAB")) {
	        		separator = '\t';
	        	}
	        	final CSVParser parser = new CSVParserBuilder().withSeparator(separator).build();
	        	if (format.contains("REACTIONCODE_CSV_WITHHEADER_")) {
	        		reader = new CSVReaderBuilder(new FileReader(fileName)).withSkipLines(1).withCSVParser(parser).build();
	        	}
	        	else {
	        		reader = new CSVReaderBuilder(new FileReader(fileName)).withCSVParser(parser).build();
	        	}
	            reader = new CSVReader(new FileReader(fileName));
	            String[] line;
	            while ((line = reader.readNext()) != null) {
	            	reactionCodes.add(line[column]);
	            }
	            reader.close();
	        } catch (IOException e) {
	            e.printStackTrace();
	        }
		}

		DecodeReactionCode decoder = new DecodeReactionCode();
		IReactionSet reactionsSet = decoder.makeReactions(reactionCodes);
		
		return reactionsSet;
	}
	
	
	/**
	 * @param fileName
	 * @param format
	 * @return
	 * @throws CDKException
	 * @throws IOException
	 * @throws CloneNotSupportedException
	 */
	public List<String> getReactionCode(String fileName, String format, boolean split) throws CDKException, IOException, CloneNotSupportedException {
		List<String> reactionCodes = new ArrayList<String>();
		
		if (format.equals("REACTIONCODE_FILE")) {
			int cpt = 0;
			reactionCount = 0;
			Path path = Paths.get(fileName);
			long lineCount = Files.lines(path).count();
			FileWriter writer = null;
			File directory = new File(new File(fileName).getParent(), "tempRCFiles/");
			if (lineCount > 10000) {
				files.clear();
				if (!directory.exists()){
			        directory.mkdir();
			    }
			}
			BufferedReader br = new BufferedReader(new FileReader(fileName)); 
			String line; 
			String reactionCode = ""; 
			while ((line = br.readLine()) != null) {
				if (line.contains(":") && line.charAt(0) == '0') {
					if (reactionCode.length() > 0) {
						if (split) {
							if (reactionCount >= 10000) {
								if (cpt == 10000) {
									if (writer != null)
										writer.close();
									File temp = new File(directory.toString() + "/chunk_" + files.size() + ".txt");
									files.add(temp);
									writer = new FileWriter(temp);
									writer.write(reactionCode + System.lineSeparator());
									cpt = 0;
									System.out.println(ColouredSystemOutPrintln.ANSI_GREEN + reactionCount + 
											" reactions are counted and a temp file has been created" + 
											ColouredSystemOutPrintln.ANSI_RESET);
								}
								else {
									writer.write(reactionCode + System.lineSeparator());
								}
							}
							else if (reactionCount < 10000) {
								reactionCodes.add(reactionCode); 
							}
							cpt++;
						}
						else {
							reactionCodes.add(reactionCode); 
						}
						reactionCount++;
					}
					reactionCode = line;
				}
				else if (line.contains(":") && line.charAt(0) != '0') {
					reactionCode += line;
				}
			}
			if (reactionCode.length() > 0) {
				cpt++;
				reactionCount++;
				if (split) {
					if (reactionCount > 10000) {
						writer.write(reactionCode + System.lineSeparator());
						if (writer != null)
							writer.close();
					}
					else {
						reactionCodes.add(reactionCode); 
					}
				}
				else {
					reactionCodes.add(reactionCode); 
				}
			}
			br.close();
			if (writer != null)
				writer.close();
		}
		else if (format.equals("REACTIONCODE_TEXT")) {
			reactionCodes.add(fileName);
		}
		else if (format.contains("REACTIONCODE_CSV_")) {
			CSVReader reader = null;
	        try {
	        	int column = Integer.parseInt(format.split("_")[3]);
	        	char separator = ','; 
	        	
	        	if (format.contains("SPACE")) {
	        		separator = ' ';
	        	}
	        	else if (format.contains("TAB")) {
	        		separator = '\t';
	        	}
	        	final CSVParser parser = new CSVParserBuilder().withSeparator(separator).build();
	        	if (format.contains("REACTIONCODE_CSV_WITHHEADER_")) {
	        		reader = new CSVReaderBuilder(new FileReader(fileName)).withSkipLines(1).withCSVParser(parser).build();
	        	}
	        	else {
	        		reader = new CSVReaderBuilder(new FileReader(fileName)).withCSVParser(parser).build();
	        	}
	            reader = new CSVReader(new FileReader(fileName));
	            String[] line;
	            while ((line = reader.readNext()) != null) {
	            	reactionCodes.add(line[column]);
	            }
	            reader.close();
	        } catch (IOException e) {
	            e.printStackTrace();
	        }
		}

		return reactionCodes;
	}
	
	public List<File> getFiles() {
		return files;
	}

	public int getReactionCount() {
		return reactionCount;
	}
	
	private void writeRDFChunk(IReactionSet reactions, String reactionFile) throws IOException {
		File directory = new File(new File(reactionFile).getParent(), "tempRDFiles/");
	    if (!directory.exists()){
	        directory.mkdir();
	    }
		File file = new File(directory.toString(), "chunk_" + files.size() + ".rdf");
		try (MDLV2000RDFWriter writer = new MDLV2000RDFWriter(new FileWriter(file))) {
			writer.setWriteAromaticBondTypes(true);
			writer.write(reactions);
			//writer.setRdFieldsProperties(map);
			writer.close();
		} catch (CDKException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		files.add(file);
	}

	public void kekulize(boolean kekulize) {
		this.kekulize = kekulize;
	}

	public void radicalize(boolean radicalize) {
		this.radicalize = radicalize;
	}
	
	
}
