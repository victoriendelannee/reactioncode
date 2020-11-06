package com.nih.main;

import static java.io.File.separator;
import static java.lang.System.out;
import static java.util.logging.Level.SEVERE;

import java.awt.Color;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.TimeoutException;
import java.util.Set;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.io.FileUtils;
import org.codehaus.jackson.map.ObjectMapper;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.Reaction;
import org.openscience.cdk.ReactionSet;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IReaction;
import org.openscience.cdk.interfaces.IReactionSet;
import org.openscience.cdk.renderer.generators.standard.StandardGenerator;
import org.openscience.cdk.tools.ILoggingTool;
import org.openscience.cdk.tools.LoggingToolFactory;

import com.nih.codes.DecodeReactionCode;
import com.nih.codes.EncodeReactionCode;
import com.nih.parser.ChemicalFormatParser;
import com.nih.reaction.PseudoMolecule;
import com.nih.tools.ColouredSystemOutPrintln;
import com.nih.tools.tools;
import com.nih.transformer.Transformer;
import com.nih.writer.MDLV2000RDFWriter;
import com.nih.writer.MDLV2000RXNWriter;
import com.opencsv.CSVWriter;

public class Main {
	
    private static final ILoggingTool LOGGER
            = LoggingToolFactory.createLoggingTool(Main.class);

	public static void main(String[] args) throws CDKException, IOException, ParseException {
		try {
			CommandLineOptions cmd = new CommandLineOptions();
			Options createEncoderOptions = cmd.createEncoderCodeOptions();
			Options createDecoderOptions = cmd.createDecoderCodeOptions();
			Options createTransformerOptions = cmd.createTransformerCodeOptions();

			DefaultParser parser1 = new DefaultParser();
			CommandLine encoderLine = parser1.parse(createEncoderOptions, args, true);
			
			DefaultParser parser2 = new DefaultParser();
			CommandLine decoderLine = parser2.parse(createDecoderOptions, args, true);
			
			DefaultParser parser3 = new DefaultParser();
			CommandLine transformerLine = parser3.parse(createTransformerOptions, args, true);

			/*
			 * Print the Header
			 */
			Helper.getHeader();

			if (encoderLine.hasOption('q')) {
				System.out.println(ColouredSystemOutPrintln.ANSI_PURPLE + "Check file format" + ColouredSystemOutPrintln.ANSI_RESET);
				ChemicalFormatParser parser = new ChemicalFormatParser();
				String file = encoderLine.getOptionValue('q');
				String format = parser.formatDetector(file);
				if (format == null) {
					System.out.println(ColouredSystemOutPrintln.ANSI_RED + "Wrong file format" + ColouredSystemOutPrintln.ANSI_RESET);
					return;
				}
				if (format.contains("REACTIONCODE")) {
					if (transformerLine.hasOption('t'))  
						transformer(transformerLine, createTransformerOptions, format);
					else
						decoder(decoderLine, createDecoderOptions, format);
				}
				else
					encoder(encoderLine, createEncoderOptions, format);
			}
			else {
				out.println("-- HELP --");
				if (encoderLine.hasOption('e')) {
					Helper.printUsageExamples();
				}
				else {
					Map<String, Options> options = new LinkedHashMap<String, Options>();
					options.put("Pseudo-Molecule & ReactionCode Encoder", createEncoderOptions);
					options.put("ReactionCode Decoder", createDecoderOptions);
					options.put("ReactionCode Transformer", createTransformerOptions);
					Helper.printHelp(options, 80, "ReactionCode", "End of Help", 5, 3, true, out);
				}
			}
		} catch (ParseException ex) {
			LOGGER.error(SEVERE, null, ex);
		} catch (Exception ex) {
			LOGGER.error(SEVERE, null, ex);
		}
	}

	private static void encoder(CommandLine commands, Options options, String format) throws CDKException, IOException, CloneNotSupportedException {
		//set default parameters for generation of reactionCode
		boolean bondType = true;
		boolean charge = true; 
		boolean hybridization = false; 
		boolean repetition = true; 
		boolean radicalize = false;
		boolean stereochemistry = true; 
		
		List<String> headers = new ArrayList<String>();
		
		//parse option and remove non-ASCII char
		String reactionFile = commands.getOptionValue("q").replaceAll("[^\\x00-\\x7F]", "");
		String outputDirectory;
		String errorsFileName;
		boolean makePseudoSmiles = false;
		boolean makeImage = false;
		boolean onlyPseudoMol = false;
		boolean writeSDF = false;
		String prefix = "encoded";
		String outputFormat = "csv";
		
		FileWriter writerSmiles = null;
		
		if (commands.hasOption('o')) {
			outputDirectory = new File(commands.getOptionValue('o')).getCanonicalPath() + separator;
		}
		else {
			if (new File(reactionFile).isFile()) {
				outputDirectory = new File(reactionFile).getParent() + separator;
			}
			else {
				outputDirectory = new File(".").getCanonicalPath() + separator;
			}
		}
		
		if (commands.hasOption('g')) {
			makePseudoSmiles = true;
			writerSmiles = new FileWriter(outputDirectory + prefix + "_pseudoSmiles.txt"); 
		}
		if (commands.hasOption('i')) {
			makeImage = true;
		}
		if (commands.hasOption('n')) {
			onlyPseudoMol = true;
		}
		if (commands.hasOption('a')) {
			bondType = false;
			hybridization = true; 
		}
		if (commands.hasOption('r')) {
			repetition = false; 
		}
		if (commands.hasOption('z')) {
			radicalize = true; 
		}
		if (commands.hasOption('s')) {
			writeSDF = true;
		}
		if (commands.hasOption('p')) {
			prefix = commands.getOptionValue('p');
		}
		if (commands.hasOption('f')) {
			outputFormat = commands.getOptionValue('f').toLowerCase();
		}
		
		errorsFileName = outputDirectory + prefix + "_ReactionCodeEncoder_errors.txt";
		removeFile(errorsFileName);
		
		//create directories
		if (writeSDF) {
			File directory = new File(outputDirectory, "SDF_PseudoMolecule");
		    if (!directory.exists()){
		        directory.mkdir();
		    }
		}
		if (makeImage) {
			File directory = new File(outputDirectory, "IMG_PseudoMolecule");
		    if (!directory.exists()){
		        directory.mkdir();
		    }
		}
		
		FileWriter outputfile = null;
		CSVWriter csvWriter = null;
		
		ObjectMapper mapper = new ObjectMapper();
		List<Object> results = new ArrayList<Object>();
		
		if (outputFormat.equals("csv")) {
			outputfile = new FileWriter(outputDirectory + prefix + "_reactionsCode.csv"); 
			csvWriter = new CSVWriter(outputfile, '\t', CSVWriter.DEFAULT_QUOTE_CHARACTER, CSVWriter.DEFAULT_ESCAPE_CHARACTER,
					CSVWriter.DEFAULT_LINE_END); 
		}
		else if (outputFormat.equals("json")) {
			outputfile = new FileWriter(outputDirectory + prefix + "_reactionsCode.json"); 
			
		}

		//parse reactions
		System.out.println(ColouredSystemOutPrintln.ANSI_BLUE + "Parse file " + ColouredSystemOutPrintln.ANSI_RESET);
		System.out.println("If the file contains more than 10000 reactions, the orginal files will be splited by chunk "
				+ "of 10000 in order to process the reactions faster");
		ChemicalFormatParser parser = new ChemicalFormatParser();
		IReactionSet reactions;

		List<File> files = new ArrayList<File>();
		
		if (format.equals("RDF_FILE")) {
			System.out.println(ColouredSystemOutPrintln.ANSI_PURPLE + "Fromat detected: RDF" 
					+ ColouredSystemOutPrintln.ANSI_RESET);
			reactions = parser.parseRDF(reactionFile, true);
			files = parser.getFiles();
			//add a null file which correspond to the reactionSet returned by parseRDF (all other reactions are 
			//written in a temp file
			files.add(0, null);
		}
		else if (format.equals("RXN_FILE")) {
			System.out.println(ColouredSystemOutPrintln.ANSI_PURPLE + "Fromat detected: RXN" 
					+ ColouredSystemOutPrintln.ANSI_RESET);
			reactions = parser.parseRXN(reactionFile);
			//add a null file which correspond to the reactionSet returned by parseRDF (all other reactions are 
			//written in a temp file
			files.add(0, null);
		}
		else if (format.equals("SMIRKS_FILE")) {
			System.out.println(ColouredSystemOutPrintln.ANSI_PURPLE + "Fromat detected: SMIRKS" 
					+ ColouredSystemOutPrintln.ANSI_RESET);
			parser.kekulize(false);
			parser.radicalize(radicalize);
			reactions = parser.parseReactionSMILES(reactionFile, 0, true);
			files = parser.getFiles();
			//add a null file which correspond to the reactionSet returned by parseRDF (all other reactions are 
			//written in a temp file
			files.add(0, null);
		}
		else if (format.equals("SMIRKS_TEXT")) {
			System.out.println(ColouredSystemOutPrintln.ANSI_PURPLE + "Fromat detected: SMIRKS" 
					+ ColouredSystemOutPrintln.ANSI_RESET);
			parser.kekulize(false);
			parser.radicalize(radicalize);
			reactions = parser.parseReactionSMILES(reactionFile);
			//add a null file which correspond to the reactionSet returned by parseRDF (all other reactions are 
			//written in a temp file
			files.add(0, null);
		}
		else {
			System.err.println("Wrong input file: file has to be formated as an RXN or an RDF file");
			return;
		}
		
		//Verify if the mapping is present for non reactionCOde input
		if (reactions.getReaction(0).getReactants().getAtomContainer(0)
				.getAtom(0).getProperty(CDKConstants.ATOM_ATOM_MAPPING) == null) {
			System.err.println("Mapping is missing in the input file");
			return;
		}
		
		int totalReactions = parser.getReactionCount();
		System.out.println(ColouredSystemOutPrintln.ANSI_PURPLE + "Number of reactions: " + totalReactions 
				+ ColouredSystemOutPrintln.ANSI_RESET);
		System.out.println(ColouredSystemOutPrintln.ANSI_BLUE + "Reactions are processing..."
				+ ColouredSystemOutPrintln.ANSI_RESET);
		
		int success = 0;
		int reactionCount = 0;
		for (int f = 0; f < files.size(); f++) {
			File file = files.get(f);
			if (files.size() > 1 && reactionCount > 0) {
				reactions = new ReactionSet();
				if (format.equals("RDF_FILE"))
					reactions = parser.parseRDF(file.toString());
				else if (format.equals("SMIRKS_FILE"))
					reactions = parser.parseReactionSMILES(file.toString(), 10000*f, false);	
			}
			for (int i = 0; i < reactions.getReactionCount(); i++) {
				IReaction reaction = reactions.getReaction(i);
				Map<Object, Object> properties = reaction.getProperties();
				IAtomContainerSet reactants = reaction.getReactants();
				IAtomContainerSet products = reaction.getProducts();
				reactionCount++;
				try {
					//aromaticity perceived when RDF is read
					for (IAtomContainer ac : reactants.atomContainers()) {
						if (format.equals("RXN_FILE") || format.equals("RDF_FILE"))
							tools.perceiveAromaticityAndflagAtomsAndBonds(ac);
						tools.attributeIDtoAtomsAndBonds(ac);
					}

					for (IAtomContainer ac : products.atomContainers()) {
						if (format.equals("RXN_FILE") || format.equals("RDF_FILE"))
							tools.perceiveAromaticityAndflagAtomsAndBonds(ac);
						tools.attributeIDtoAtomsAndBonds(ac);
					}
				}
				catch (Exception e) {
					dataFile(reaction.getID() + "\t The mapping information is probably wrong. Please check that each atom "
							+ "has value superior to 0.\n", errorsFileName);
					continue;
				}
				
				//Annotate the reaction with right constants
				PseudoMolecule pm = new PseudoMolecule();
				IAtomContainer pseudoMolecule = null;

				try {
					pm.reactionAnnotator(reaction);
					pseudoMolecule = pm.makePseudoMolecule(reactants, products);
					pseudoMolecule.setProperty(CDKConstants.TITLE, "PSEUDO-MOLECULE");
				}
				catch (Exception e) {
					dataFile(reaction.getID() + "\t Something get wrong during the PseudoMoleclule generation\n", errorsFileName);
					continue;
				}

				if (writeSDF == true) {
					try {
						pm.writePseudoMoleculeSDFile(prefix+"_"+reactionCount, new File(outputDirectory, "SDF_PseudoMolecule").toString(), 
								pseudoMolecule.clone());
					}
					catch (Exception e) {
						dataFile(reaction.getID() + "\t The SDF could not be generated\n", errorsFileName);
					}
				}
				if (makeImage == true) {
					try {
						pm.writePseudoMoleculeImage(new File(outputDirectory, "IMG_PseudoMolecule").toString(), prefix+"_"+reactionCount);
					}
					catch (Exception e) {
						dataFile(reaction.getID() + "\t The IMG could not be generated\n", errorsFileName);
					}
				}
				
				String pseudoSmiles = null;

				try {
					pseudoSmiles = pm.GetPseudoMoleculeSmiles(pseudoMolecule.clone());
				}
				catch (Exception e) {
					dataFile(reaction.getID() + "\t The pseudoSmiles could not be generated\n", errorsFileName);
				}
								
				if (makePseudoSmiles) {
					if (pseudoSmiles != null) {
						writerSmiles.write(pseudoSmiles+"\n");
					}
				}

				if (onlyPseudoMol == false) {
					//get atom duplication number
					Map<String, Integer> numberOfRepetitions = pm.atomRepetition();
					Set<IAtom> reactionCenterAtom = pm.getReactioncenter();

					if (reactionCenterAtom.size() == 0) {
						dataFile(reaction.getID() + "\t No atom has been identified as a reaction centre. The mapping is probably incorrect\n", errorsFileName);
						continue;
					}
					try {
						//make reactionCode
						//set parameters for generation of reactionCode
						EncodeReactionCode encoder = new EncodeReactionCode();
						
						//set parameters for generation of reactionCode
						encoder.setBondType(bondType);
						encoder.setCharge(charge);
						encoder.setHybridization(hybridization);
						encoder.setRepetition(repetition);
						encoder.setStereochemistry(stereochemistry);
						Map<String,String> reactionCodeMap = encoder.generateReactionCode(reactionCenterAtom, reactants, 
								pseudoMolecule, numberOfRepetitions);
						//String reactionCode = ReactionCode.reactionCodeMapToStringMultiLines(reactionCodeMap);
						String reactionCode = encoder.reactionCodeMapToStringOneLine(reactionCodeMap);
	
						if (outputFormat.equals("json")) {
							//set JSON object
							Map<String,String> info = new HashMap<String,String>();
							info.put("id", reaction.getID());
							info.put("pseudosmiles", pseudoSmiles);
							for (Entry<Object, Object> e : properties.entrySet()) {
								info.put(e.getKey().toString(), e.getValue().toString());
							}
							info.put("reactionCode", reactionCode);
							results.add(info);
						}
						else if (outputFormat.equals("csv")) {
							try { 
								if (success == 0) {
									String[] header = new String[properties.size()+3];
									String[] data = new String[properties.size()+3];
									header[0] = "id";
									headers.add("id");
									data[0] = reaction.getID();
									header[1] = "pseudosmiles";
									headers.add("pseudosmiles");
									data[1] = pseudoSmiles;
									int cpt = 2;
									for (Entry<Object, Object> e : properties.entrySet()) {
										String propertyName = e.getKey().toString();
										header[cpt] = propertyName;
										headers.add(propertyName);
										Object value = e.getValue();
										if (value == null) {
											data[cpt] = "NA";
										}
										else {
											data[cpt] = value.toString();
										}
										cpt++;
									}
									header[cpt] = "reactionCode";
									headers.add("reactionCode");
									data[cpt] = reactionCode;
									csvWriter.writeNext(header);
									csvWriter.writeNext(data);
								}
								else {
									String[] data = new String[properties.size()+3];
									data[0] = reaction.getID();
									data[1] = pseudoSmiles;
									int cpt = 2;
									for (Entry<Object, Object> e : properties.entrySet()) {
										String propertyName = e.getKey().toString();
										if (headers.contains(propertyName)) {
											Object value = e.getValue();
											if (value == null) {
												data[cpt] = "NA";
											}
											else {
												data[cpt] = value.toString();
											}
											cpt++;
										}
									}
									data[cpt] = reactionCode;
									csvWriter.writeNext(data);
								}
							}
							catch (Exception e) { 
								// TODO Auto-generated catch block 
								e.printStackTrace(); 
							}
						}
						success++;
					}
					catch (Exception e) {
						dataFile(reaction.getID() + "\t Something went wrong during the ReactionCode generation\n", errorsFileName);
					}
				}
				if (reactionCount % 1000 == 0) {
					System.out.println(ColouredSystemOutPrintln.ANSI_GREEN + reactionCount + "/" +
							totalReactions + " reactions are processed..." + ColouredSystemOutPrintln.ANSI_RESET);
				}
			}
		}

		System.out.println(ColouredSystemOutPrintln.ANSI_PURPLE + "All reactions are processed...");

		System.out.println(ColouredSystemOutPrintln.ANSI_BLUE + "Output files are finishing to write...");
		if (outputFormat.equals("json")) {
			//write JSON
			try {
				mapper.writerWithDefaultPrettyPrinter().writeValue(outputfile, results);
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
		else if (outputFormat.equals("csv")) {
			csvWriter.close();
		}
		
		if (makePseudoSmiles) {
			writerSmiles.close();
		}
		
		if (new File(outputDirectory, "tempRDFiles").exists()) {
			FileUtils.deleteDirectory(new File(outputDirectory, "tempRDFiles"));
		}
		
		int fails = totalReactions-success;
		System.out.println(ColouredSystemOutPrintln.ANSI_GREEN + "Done ;)");
		System.out.println(ColouredSystemOutPrintln.ANSI_GREEN +
				success + " reactions have been process with success");
		System.out.println(ColouredSystemOutPrintln.ANSI_RED +
				fails + " reactions have failed"+ ColouredSystemOutPrintln.ANSI_RESET);

	}
	
	private static void decoder(CommandLine commands, Options options, String format) throws CDKException, IOException, CloneNotSupportedException {

		//parse option
		String reactionFile = commands.getOptionValue("q").replaceAll("[^\\x00-\\x7F]", "");
		String outputDirectory;
		String errorsFileName;
		boolean makeImage = false;
		boolean addHydrogenForLastLayer = false;
		boolean writeSMIRKS = false;
		boolean writeSMIRKSWMapping = false;
		boolean writeRXN = false;
		boolean writeRDF = false;
		String prefix = "decoded";

		if (commands.hasOption('o')) {
			outputDirectory = new File(commands.getOptionValue('o')).getCanonicalPath() + separator;
		}
		else {
			if (new File(reactionFile).isFile()) {
				outputDirectory = new File(reactionFile).getParent() + separator;
			}
			else {
				outputDirectory = new File(".").getCanonicalPath() + separator;
			}
		}

		if (commands.hasOption('i')) {
			makeImage = true;
		}
		if (commands.hasOption('h')) {
			addHydrogenForLastLayer = false;
		}
		if (commands.hasOption('s')) {
			writeSMIRKS = true;
		}
		if (commands.hasOption('m')) {
			writeSMIRKSWMapping = true;
		}
		if (commands.hasOption('r')) {
			writeRDF = true;
		}
		if (commands.hasOption('x')) {
			writeRXN = true;
		}
		if (commands.hasOption('p')) {
			prefix = commands.getOptionValue('p');
		}
		
		if (commands.getOptions().length < 2 || (commands.getOptions().length == 2 && commands.hasOption('p'))) {
			System.err.println(ColouredSystemOutPrintln.ANSI_RED +
					"The command line has to have the option -q and at least one of these options: -s, -m, -i, -r, -x"
					+ ColouredSystemOutPrintln.ANSI_RESET);
			return;
		}
		
		errorsFileName = outputDirectory + prefix + "_ReactionCodeDecoder_errors.txt";
		removeFile(errorsFileName);

		//create directories
		if (writeRXN) {
			File directory = new File(outputDirectory, "RXN_Reaction");
			if (!directory.exists()){
				directory.mkdir();
			}
		}
		if (makeImage) {
			File directory = new File(outputDirectory, "IMG_ReactionDecoded");
			if (!directory.exists()){
				directory.mkdir();
			}
		}

		//delete file with same name
		if (new File(outputDirectory + prefix + "_reactionSMILES.txt").exists()) {
			removeFile(outputDirectory + prefix + "_reactionSMILES.txt");
		}
		if (new File(outputDirectory + prefix + "_SMIRKS.txt").exists()) {
			removeFile(outputDirectory + prefix + "_SMIRKS.txt");
		}

		//parse reactions
		System.out.println("Decode ReactionCode");
		ChemicalFormatParser parser = new ChemicalFormatParser();
		List<String> reactionCodes = parser.getReactionCode(reactionFile, format, true);
		
		List<File> files = new ArrayList<File>();
		files = parser.getFiles();
		//add a null file which correspond to the list returned by getReactionCode (all other code are 
		//written in a temp file
		files.add(0, null);

		System.out.println(ColouredSystemOutPrintln.ANSI_GREEN +
				parser.getReactionCount() + " reactions are decoding "+ ColouredSystemOutPrintln.ANSI_RESET);

		DecodeReactionCode decoder = new DecodeReactionCode();
		decoder.setAddHydrogenForLastLayer(addHydrogenForLastLayer); 
		IReactionSet reactions = DefaultChemObjectBuilder.getInstance().newInstance(IReactionSet.class);
		int cpt = 1;
		int reactionCount = 0;
		int filecpt = 1;
		int success = 0;
		int fails = 0;
		boolean failed = false;
		int id = 0;
		for (int i = 0; i < files.size(); i++) {
			File tfile = files.get(i);
			if (files.size() > 1 && reactionCount > 0) {
				reactionCodes.clear();
				reactionCodes = parser.getReactionCode(tfile.toString(), format, false);
			}
			for (String reactionCode : reactionCodes) {
				IReaction reaction;
				reactionCount++;
				try {
					reaction = decoder.decode(reactionCode);
					reaction.setID(id+"");
				}
				catch (Exception e){
					reaction = new Reaction();
					reaction.setID(id+"");
					dataFile(id + "\t decoding failure \t" + reactionCode + "\n", errorsFileName);
					failed = true;
				} 

				if (reaction.getProperty("mappingError") != null) {
					//write failed reactions
					dataFile(id + "\t mapping error \t" + reactionCode + "\n", errorsFileName);
					failed = true;
				}
				reactions.addReaction(reaction);
				id++;
				cpt++;

				if (writeRXN || makeImage || writeSMIRKS || writeSMIRKSWMapping) {
					if (writeRXN) {
						File file = new File(outputDirectory + "/RXN_Reaction/reaction_" + prefix + " " + id + "_decodedReactions.rxn");
						try (MDLV2000RXNWriter writer = new MDLV2000RXNWriter(new FileWriter(file))) {
							writer.setWriteAromaticBondTypes(true);
							writer.write(reactions);
							//writer.setRdFieldsProperties(map);
							writer.close();
						} catch (CDKException e) {
							// TODO Auto-generated catch block
							e.printStackTrace();
							dataFile(id + "\t Can not write a proper RXN \t" + reactionCode + "\n", errorsFileName);
						}
					}
					if (writeSMIRKS){
						try {
							dataFile(id + "\t" + tools.makeSmiles(reaction, false) + "\n", outputDirectory + prefix + "_reactionSMILES.txt");
						}
						catch (Exception e) {
							e.printStackTrace();
						}
					}
					if (writeSMIRKSWMapping) {
						try {
							dataFile(id + "\t" + tools.makeSmiles(reaction, true) + "\n", outputDirectory + prefix + "_SMIRKS.txt");
						}
						catch (Exception e) {
							e.printStackTrace();
							dataFile(id + "\t Can not write a proper SMIRKS \t" + reactionCode + "\n", errorsFileName);
						}
					}
					if (makeImage) {
						try {
							//Convert IQueryAtomContainer to AtomContainer
							IReaction newReaction = tools.convertReactionContainingIQueryAtomContainer(reaction);
							//force delocalised in order to mark aromatic as dash bond (even if the ring is not complete)
							new org.openscience.cdk.depict.DepictionGenerator()
							.withHighlight(newReaction.getProperty("bondsFormed"), Color.BLUE)
							.withHighlight(newReaction.getProperty("bondsCleaved"), Color.RED)
							.withHighlight(newReaction.getProperty("bondsOrder"), Color.GREEN)
							.withHighlight(newReaction.getProperty("reactionCenter"), Color.LIGHT_GRAY)
							.withOuterGlowHighlight().withAtomColors().withAtomMapNumbers()
							.depict(newReaction).writeTo(outputDirectory + "/IMG_ReactionDecoded/reaction_" + reaction.getID() + ".pdf");
						}
						catch (Exception e) {
							e.printStackTrace();
							dataFile(id + "\t Can not depict the reaction \t" + reactionCode + "\n", errorsFileName);
						}

					}
				}
				
				if (cpt == 5000) {
					System.out.println(ColouredSystemOutPrintln.ANSI_GREEN + 5000*filecpt +
							" have been decoded "+ ColouredSystemOutPrintln.ANSI_RESET);
					if (writeRDF) {
						System.out.println("RD file is being wrting...");
						File file = new File(outputDirectory + filecpt + "_decodedReactions.rdf");
						try (MDLV2000RDFWriter writer = new MDLV2000RDFWriter(new FileWriter(file))) {
							writer.setWriteAromaticBondTypes(true);
							writer.write(reactions);
							//writer.setRdFieldsProperties(map);
							writer.close();
						} catch (CDKException e) {
							// TODO Auto-generated catch block
							e.printStackTrace();
							dataFile(id + "\t Can not write a proper RDF \t" + reactionCode + "\n", errorsFileName);
						}
						System.out.println(ColouredSystemOutPrintln.ANSI_GREEN +
								"RD file writing is complete "+ ColouredSystemOutPrintln.ANSI_RESET);
					}
					reactions.removeAllReactions();
					cpt = 1;
					filecpt++;
				}

				if (failed) 
					fails += 1;
				else
					success += 1;
				failed = false;

			}
		}
		
		if (new File(outputDirectory, "tempRCFiles").exists()) {
			FileUtils.deleteDirectory(new File(outputDirectory, "tempRCFiles"));
		}
		
		System.out.println(ColouredSystemOutPrintln.ANSI_GREEN +
				success + " reactions have been process with success");
		System.out.println(ColouredSystemOutPrintln.ANSI_RED +
				fails + " reactions have failed"+ ColouredSystemOutPrintln.ANSI_RESET);
		System.out.println(ColouredSystemOutPrintln.ANSI_GREEN + "Done ;)" + ColouredSystemOutPrintln.ANSI_RESET);
	}

	private static void transformer(CommandLine commands, Options options, String format) throws CDKException, IOException, CloneNotSupportedException {
		//parse option
		String reactionFile = commands.getOptionValue('q').replaceAll("[^\\x00-\\x7F]", "");
		String reagents = commands.getOptionValue('t').replaceAll("[^\\x00-\\x7F]", "");
		String outputDirectory;
		String errorsFileName;
		boolean makeImage = false;
		boolean writeSMIRKS = false;
		boolean writeSMIRKSWMapping = false;
		boolean writeRXN = false;
		boolean writeRDF = false;
		boolean multipleReactions = false;
		int stoichiometry = 1;
		boolean checkValence = false;
		int time = 60;
		String prefix = "transformed";
		
		if (commands.hasOption('o')) {
			outputDirectory = new File(commands.getOptionValue('o')).getCanonicalPath() + separator;
		}
		else {
			if (new File(reactionFile).isFile()) {
				outputDirectory = new File(reactionFile).getParent() + separator;
			}
			else {
				outputDirectory = new File(".").getCanonicalPath() + separator;
			}
		}
		if (commands.hasOption('a')) {
			multipleReactions = true;
		}
		if (commands.hasOption('z')) {
			stoichiometry = Integer.parseInt(commands.getOptionValue('z'));
		}
		if (commands.hasOption('z')) {
			time = Integer.parseInt(commands.getOptionValue('z'));
		}
		if (commands.hasOption('c')) {
			checkValence = true;
		}
		if (commands.hasOption('i')) {
			makeImage = true;
		}
		if (commands.hasOption('s')) {
			writeSMIRKS = true;
		}
		if (commands.hasOption('m')) {
			writeSMIRKSWMapping = true;
		}
		if (commands.hasOption('r')) {
			writeRDF = true;
		}
		if (commands.hasOption('x')) {
			writeRXN = true;
		}
		if (commands.hasOption('p')) {
			prefix = commands.getOptionValue('p');
		}
		
		if (commands.getOptions().length < 2 || (commands.getOptions().length == 2 && commands.hasOption('p'))) {
			System.err.println(ColouredSystemOutPrintln.ANSI_RED +
					"The command line has to have the option -q and at least one of these options: -s, -m, -i, -r, -x"
					+ ColouredSystemOutPrintln.ANSI_RESET);
			return;
		}
		
		errorsFileName = outputDirectory + prefix + "_ReactionCodeTransformer_errors.txt";
		removeFile(errorsFileName);
		
		

		//create directories
		if (writeRXN) {
			File directory = new File(outputDirectory, "RXN_ReactionTransformed");
			if (!directory.exists()){
				directory.mkdir();
			}
		}
		if (makeImage) {
			File directory = new File(outputDirectory, "IMG_ReactionTransformed");
			if (!directory.exists()){
				directory.mkdir();
			}
		}

		//delete file with same name
		if (new File(outputDirectory + prefix + "_reactionSMILES.txt").exists()) {
			removeFile(outputDirectory + prefix + "_reactionSMILES.txt");
		}
		if (new File(outputDirectory + prefix + "_SMIRKS.txt").exists()) {
			removeFile(outputDirectory + prefix + "_SMIRKS.txt");
		}
		if (new File(outputDirectory + prefix + "_reactionSMILES_1.txt").exists()) {
			removeFile(outputDirectory + prefix + "_reactionSMILES_1.txt");
		}
		if (new File(outputDirectory + prefix + "_SMIRKS_1.txt").exists()) {
			removeFile(outputDirectory + prefix + "_SMIRKS_1.txt");
		}
		
		//parse reactions
		System.out.println("Parse transformation");
		ChemicalFormatParser parser = new ChemicalFormatParser();
		List<String> reactionCodes = parser.getReactionCode(reactionFile, format, false);
		
		//List<File> files = new ArrayList<File>();
		//files = parser.getFiles();
		//add a null file which correspond to the list returned by getReactionCode (all other code are 
		//written in a temp file
		//files.add(0, null);

		System.out.println(ColouredSystemOutPrintln.ANSI_GREEN +
				reactionCodes.size() + " transformations are beeing applied "+ ColouredSystemOutPrintln.ANSI_RESET);

		//check if it is SMILES/SMART or a file containing path of SDF or SMILES/SMART
		String reagentsFormat = "SMILES";
		File reagentsFile = new File(reagents); 
		List<String> reagentsList = new ArrayList<String>(); 
		if (reagentsFile.isFile()) {
			BufferedReader br = new BufferedReader(new FileReader(reagentsFile)); 
			String line; 
			boolean firstLine = true;
			while ((line = br.readLine()) != null) {
				if (firstLine) {
					try {
						parser.parseSMILES(line, false);
						reagentsFormat = "setOfSMILES";
					}
					catch (Exception e) {
						if (new File(line).isFile()) {
							reagentsFormat = "setOfSDF";
						}
					}
					firstLine = false;
				}
				reagentsList.add(line);
			}
		}

		Transformer transformer = new Transformer();
		transformer.setStoichiometry(stoichiometry);
		transformer.setCheckValence(checkValence);
		
		DecodeReactionCode decoder = new DecodeReactionCode();
		decoder.setCalculateExpr(true);
		decoder.setAddHydrogenForLastLayer(false);
		IReactionSet reactions = DefaultChemObjectBuilder.getInstance().newInstance(IReactionSet.class);
		int cpt = 1;
		int reactionCount = 0;
		int filecpt = 1;
		int success = 0;
		int fails = 0;
		int id = 0;
		/*
		for (int i = 0; i < files.size(); i++) {
			File tfile = files.get(i);
			if (files.size() > 1 && reactionCount > 0) {
				reactionCodes.clear();
				reactionCodes = parser.getReactionCode(tfile.toString(), format);
			}
			*/
		List<IAtomContainerSet> allReactants = new ArrayList<IAtomContainerSet>();
		if (reagentsFile.isFile()) {
			BufferedReader br = new BufferedReader(new FileReader(reagentsFile)); 
			String line; 
			while ((line = br.readLine()) != null) {
				IAtomContainerSet reactants = null;
				if (reagentsFormat.equals("setOfSMILES")){
					reactants = parser.parseSMILES(line, false);
					allReactants.add(reactants);
				}
				else if (reagentsFormat.equals("setOfSMILES")){
					reactants = parser.parseSDF(line).getReactants();
					allReactants.add(reactants);
				}
			}
		}
			for (int i = 0; i< reactionCodes.size(); i++) {
				String reactionCode = reactionCodes.get(i);
				reactionCount++;
				IReaction transform;
				try {
					transform = decoder.decode(reactionCode);
					transform.setID(id+"");
				}
				catch (Exception e){
					transform = new Reaction();
					transform.setID(id+"");
					dataFile(id + "\t decoding failure \t" + reactionCode + "\n", errorsFileName);
				} 
	
				if (transform.getProperty("mappingError") != null) {
					//write failed reactions
					dataFile(id + "\t mapping error \t" + reactionCode + "\n", errorsFileName);
				}
	//			File reagentsFile = new File(reagents); 
				List<IReaction> generatedReactions = new ArrayList<IReaction>();
				IReaction transformCopy = transform;
			if (reagentsFormat.equals("SMILES")) {
				final ExecutorService service = Executors.newSingleThreadExecutor();
				try {
		            final Future<List<IReaction> > f = service.submit(() -> {
		                return transformer.transform(reagents, transformCopy);
		            });

		            generatedReactions = f.get(time, TimeUnit.SECONDS);
		        } catch (final TimeoutException e) {
		            System.err.println("Calculation took to long");
		        } catch (final Exception e) {
		            throw new RuntimeException(e);
		        } finally {
		        	cpt++;
		        	id++;
		            service.shutdown();
		        }
			} 
			else if (reagentsFormat.equals("setOfSMILES")) {
				IAtomContainerSet reactants = allReactants.get(i);
				final ExecutorService service = Executors.newSingleThreadExecutor();
				try {
		            final Future<List<IReaction> > f = service.submit(() -> {
		                return transformer.transform(reactants, transformCopy);
		            });

		            generatedReactions = f.get(time, TimeUnit.SECONDS);
		        } catch (final TimeoutException e) {
		            System.err.println("Calculation took to long");
		        } catch (final Exception e) {
		            throw new RuntimeException(e);
		        } finally {
		        	cpt++;
		        	id++;
		            service.shutdown();
		        }
			} 
			else if (reagentsFormat.equals("setOfSDF")) {
				IAtomContainerSet reactants = allReactants.get(i);
				final ExecutorService service = Executors.newSingleThreadExecutor();
				try {
		            final Future<List<IReaction> > f = service.submit(() -> {
		                return transformer.transform(reactants, transformCopy);
		            });

		            generatedReactions = f.get(time, TimeUnit.SECONDS);
		        } catch (final TimeoutException e) {
		            System.err.println("Calculation took to long");
		        } catch (final Exception e) {
		            throw new RuntimeException(e);
		        } finally {
		        	cpt++;
		        	id++;
		            service.shutdown();
		        }
			} 
			else {
				System.err.println("invalid format of Reagents or does not match with the reactionCodes number");
				return;
			}
			for (IReaction r : generatedReactions) {
				r.setID(transform.getID());
				reactions.addReaction(r);
			}
			
			if (generatedReactions.size() == 0) 
				fails += 1;
			else
				success += 1;
				
				if (writeRXN || makeImage || writeSMIRKS || writeSMIRKSWMapping) {
					if (writeRXN) {
						for (int j = 0; j < generatedReactions.size(); j++) {
							IReaction reaction = generatedReactions.get(j);
							File file;
							if (multipleReactions){
								file = new File(outputDirectory + "/RXN_ReactionTransformed/reaction_" + prefix + " " + id + "_" + j + "_transformedReactions.rxn");
							}
							else{
								file = new File(outputDirectory + "/RXN_ReactionTransformed/reaction_" + prefix + " " + id + "_transformedReactions.rxn");
							}
	
							try (MDLV2000RXNWriter writer = new MDLV2000RXNWriter(new FileWriter(file))) {
								writer.setWriteAromaticBondTypes(true);
								writer.write(reaction);
								//writer.setRdFieldsProperties(map);
								writer.close();
							} catch (CDKException e) {
								// TODO Auto-generated catch block
								e.printStackTrace();
							}
						}
					}
					if (writeSMIRKS){
						String smiles = "";
						for (IReaction reaction : generatedReactions) {
							smiles += tools.makeSmiles(reaction, false) + "\t";
						}
						if (generatedReactions.isEmpty()) {
							smiles = "no or too many possible products";
						}
						try {
							dataFile(id-1 + "\t" + smiles + "\n", outputDirectory + prefix + "_reactionSMILES_" + filecpt + ".txt");
						}
						catch (Exception e) {
							e.printStackTrace();
						}
					}
					if (writeSMIRKSWMapping) {
						String smiles = "";
						for (IReaction reaction : generatedReactions) {
							smiles += tools.makeSmiles(reaction, true) + "\t";
						}
						if (generatedReactions.isEmpty()) {
							smiles = "no or too many possible products";
						}
						try {
							dataFile(id-1 + "\t" + smiles + "\n", outputDirectory + prefix + "_SMIRKS_" + filecpt + ".txt");
						}
						catch (Exception e) {
							e.printStackTrace();
						}
					}
					if (makeImage) {
						try {
							//force delocalised in order to mark aromatic as dash bond (even if the ring is not complete)
							for (int j = 0; j < generatedReactions.size(); j++) {
								IReaction reaction = generatedReactions.get(j);
								new org.openscience.cdk.depict.DepictionGenerator().withParam(StandardGenerator.ForceDelocalisedBondDisplay.class, true)
								//.withHighlight(reaction.getProperty("bondsFormed"), Color.BLUE)
								//.withHighlight(reaction.getProperty("bondsCleaved"), Color.RED)
								//.withHighlight(reaction.getProperty("bondsOrder"), Color.GREEN)
								//.withHighlight(reaction.getProperty("reactionCenter"), Color.LIGHT_GRAY)
								.withOuterGlowHighlight().withAtomColors().withAtomMapNumbers()
								.depict(reaction).writeTo(outputDirectory + "/IMG_ReactionTransformed/reaction_" + reaction.getID() + "_" + j + ".pdf");
							}
						}
						catch (Exception e) {
							e.printStackTrace();
						}
	
					}
				}

				if (cpt >= 250) {
					System.out.println(ColouredSystemOutPrintln.ANSI_GREEN + 250*filecpt +
							" have been transformed "+ ColouredSystemOutPrintln.ANSI_RESET);
					if (writeRDF) {
						System.out.println("RD file is being writing...");
						File file = new File(outputDirectory + filecpt + "_transformedReactions.rdf");
						try (MDLV2000RDFWriter writer = new MDLV2000RDFWriter(new FileWriter(file))) {
							writer.setWriteAromaticBondTypes(true);
							writer.write(reactions);
							//writer.setRdFieldsProperties(map);
							writer.close();
						} catch (CDKException e) {
							// TODO Auto-generated catch block
							e.printStackTrace();
						}
						System.out.println(ColouredSystemOutPrintln.ANSI_GREEN +
								"RD file writing is complete "+ ColouredSystemOutPrintln.ANSI_RESET);
					}
					reactions.removeAllReactions();
					cpt = 1;
					filecpt++;
					if (new File(outputDirectory + prefix + "_reactionSMILES_" + filecpt + ".txt").exists()) {
						removeFile(outputDirectory + prefix + "_reactionSMILES_" + filecpt + ".txt");
					}
					if (new File(outputDirectory + prefix + "_SMIRKS_" + filecpt + ".txt").exists()) {
						removeFile(outputDirectory + prefix + "_SMIRKS_" + filecpt + ".txt");
					}
				}
			}
		//}
		System.out.println(ColouredSystemOutPrintln.ANSI_GREEN +
				success + " reactions have been process with success");
		System.out.println(ColouredSystemOutPrintln.ANSI_RED +
				fails + " reactions have failed"+ ColouredSystemOutPrintln.ANSI_RESET);
		System.out.println(ColouredSystemOutPrintln.ANSI_GREEN + "Done ;)" + ColouredSystemOutPrintln.ANSI_RESET);
	}

	
	/*private static void decoder(CommandLine commands, Options options, String format) throws CDKException, IOException, CloneNotSupportedException {
		
		//parse option
		String reactionFile = commands.getOptionValue("q").replaceAll("[^\\x00-\\x7F]", "");
		String outputDirectory;
		String errorsFileName;
		boolean makeImage = false;
		boolean writeSMIRKS = false;
		boolean writeSMIRKSWMapping = true;
		boolean writeRXN = false;
		boolean writeRDF = false;
		String prefix = "PR";
		
		if (commands.hasOption('o')) {
			outputDirectory = new File(commands.getOptionValue('o')).getCanonicalPath() + separator;
		}
		else {
			if (new File(reactionFile).isFile()) {
				outputDirectory = new File(reactionFile).getParent() + separator;
			}
			else {
				outputDirectory = new File(".").getCanonicalPath() + separator;
			}
		}
		
		errorsFileName = outputDirectory + "ReactionCodeDecoder_errors.txt";
		removeFile(errorsFileName);
		
		if (commands.hasOption('i')) {
			makeImage = true;
		}
		if (commands.hasOption('s')) {
			writeSMIRKS = true;
		}
		if (commands.hasOption('s')) {
			writeSMIRKSWMapping = true;
		}
		if (commands.hasOption('r')) {
			writeRDF = true;
		}
		if (commands.hasOption('x')) {
			writeRXN = true;
		}
		if (commands.hasOption('p')) {
			prefix = commands.getOptionValue('p');
		}

		//create directories
		if (writeRXN) {
			File directory = new File(outputDirectory, "RXN_Reaction");
		    if (!directory.exists()){
		        directory.mkdir();
		    }
		}
		if (makeImage) {
			File directory = new File(outputDirectory, "IMG_ReactionDecoded");
		    if (!directory.exists()){
		        directory.mkdir();
		    }
		}

		//parse reactions
		System.out.println("Decode ReactionCode");
		ChemicalFormatParser parser = new ChemicalFormatParser();
		IReactionSet reactions = parser.parseReactionCode(reactionFile, format);
		
		System.out.println(ColouredSystemOutPrintln.ANSI_GREEN +
				reactions.getReactionCount() + " reactions have been decoded"+ ColouredSystemOutPrintln.ANSI_RESET);

		if (writeRDF) {
			System.out.println("RD file is being writing...");
			File file = new File(outputDirectory + "decodedReactions.rdf");
			try (MDLV2000RDFWriter writer = new MDLV2000RDFWriter(new FileWriter(file))) {
				writer.write(reactions);
				//writer.setRdFieldsProperties(map);
				writer.close();
			} catch (CDKException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			System.out.println(ColouredSystemOutPrintln.ANSI_GREEN +
					"RD file writing is complete "+ ColouredSystemOutPrintln.ANSI_RESET);
		}
		if (writeRXN || makeImage || writeSMIRKS || writeSMIRKSWMapping) {
			FileWriter smirksWriter = null;
			FileWriter smirksWriterWithMapping = null;
			if (writeSMIRKS){
				smirksWriter = new FileWriter(outputDirectory + prefix + "_reactionSMILES.txt"); 
			}
			if (writeSMIRKSWMapping) {
				smirksWriterWithMapping = new FileWriter(outputDirectory + prefix + "_reactionSMILESAndMapping.txt"); 
			}
			System.out.println("Files (RXN, Reaction SMILES, and Images) are beeing writing...");
			for (int i = 0; i < reactions.getReactionCount(); i++) {
				IReaction reaction = reactions.getReaction(i);
				if (writeRXN) {
					File file = new File(outputDirectory + "/RXN_Reaction/", "reaction_" + i + ".rxn");
					try (MDLV2000RXNWriter writer = new MDLV2000RXNWriter(new FileWriter(file))) {
						writer.write(reaction);
						//writer.setRdFieldsProperties(map);
						writer.close();
					} catch (CDKException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
				}
				if (makeImage) {
					try {
						new org.openscience.cdk.depict.DepictionGenerator()
						.withHighlight(reaction.getProperty("bondsFormed"), Color.BLUE)
						.withHighlight(reaction.getProperty("bondsCleaved"), Color.RED)
						.withHighlight(reaction.getProperty("bondsOrder"), Color.GREEN)
						.withHighlight(reaction.getProperty("reactionCenter"), Color.LIGHT_GRAY)
						.withOuterGlowHighlight().withAtomColors().withAtomMapNumbers()
						.depict(reaction).writeTo(outputDirectory + "/IMG_ReactionDecoded/reaction_" + i + ".pdf");
					}
					catch (Exception e) {
						e.printStackTrace();
					}
					
				}
				if (writeSMIRKSWMapping) {
					try {
						smirksWriter.write(tools.makeSmiles(reaction, false) + "\n");
					}
					catch (Exception e) {
						e.printStackTrace();
					}
					smirksWriterWithMapping.write(tools.makeSmiles(reaction, true) + "\n");
				}
				if (writeSMIRKS) {
					try {
						smirksWriter.write(tools.makeSmiles(reaction, false) + "\n");
					}
					catch (Exception e) {
						e.printStackTrace();
					}
				}
			}
			
			//write failed reactions
			dataFile((List<String>)reactions.getProperty("errors"), errorsFileName);
			
			System.out.println(ColouredSystemOutPrintln.ANSI_GREEN +
					"Files writing is complete "+ ColouredSystemOutPrintln.ANSI_RESET);
			if (writeSMIRKS)
				smirksWriter.close(); 
			if (writeSMIRKSWMapping)
				smirksWriterWithMapping.close(); 
		}
		System.out.println(ColouredSystemOutPrintln.ANSI_GREEN + "Done ;)" + ColouredSystemOutPrintln.ANSI_RESET);
	}*/

	
	private static void dataFile(String data, String filename) {
		BufferedWriter bw = null;
		FileWriter fw = null;
		try {
			File file = new File(filename);
			// if file doesnt exists, then create it
			if (!file.exists()) {
				file.createNewFile();
			}
			// true = append file
			fw = new FileWriter(file.getAbsoluteFile(), true);
			bw = new BufferedWriter(fw);
			bw.write(data);
			//System.out.println("Done");
		} catch (IOException e) {
			e.printStackTrace();
		} finally {
			try {
				if (bw != null)
					bw.close();
				if (fw != null)
					fw.close();
			} catch (IOException ex) {
				ex.printStackTrace();
			}
		}
	}
	
	private static void dataFile(List<String> data, String filename) {
		BufferedWriter bw = null;
		FileWriter fw = null;
		try {
			File file = new File(filename);
			// if file doesnt exists, then create it
			if (!file.exists()) {
				file.createNewFile();
			}
			// true = append file
			fw = new FileWriter(file.getAbsoluteFile(), true);
			bw = new BufferedWriter(fw);
			for (String line : data) {
				bw.write(line);
			}
			//System.out.println("Done");
		} catch (IOException e) {
			e.printStackTrace();
		} finally {
			try {
				if (bw != null)
					bw.close();
				if (fw != null)
					fw.close();
			} catch (IOException ex) {
				ex.printStackTrace();
			}
		}
	}
	
	private static void removeFile(String filename) {
		try {
			File file = new File(filename);
			if (file.exists()) {
				file.delete();
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
}
