package com.nih.main;

import static java.io.File.separator;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.Map.Entry;

import org.codehaus.jackson.map.ObjectMapper;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IReaction;
import org.openscience.cdk.interfaces.IReactionSet;

import com.nih.codes.EncodeReactionCode;
import com.nih.parser.ChemicalFormatParser;
import com.nih.reaction.PseudoMolecule;
import com.nih.tools.tools;
import com.opencsv.CSVWriter;

public class TestEncoder {

	public static void main(String[] args) throws CDKException, IOException {
		// TODO Auto-generated method stub

		//"/datahdd/SAVI/Alexey/testMapping/testMapping2.rdf"

		//parse option
		//String reactionFile = "/datahdd/SAVI/Alexey/testMapping/testMapping3.rdf";
		//String reactionFile = "/home/delanneev/Downloads/devendraFromated.txt";
		//String reactionFile = "[CH3:1][Si:4]([Cl:3])([CH3:6])[CH3:8].[CH2:2]([C:5]=1[CH:11]=[C:13]([NH2:15])[CH:16]=[CH:14][CH:12]1)[C:7](=[O:9])[OH:10]>>[CH2:1]([CH:2]([C:5]=1[CH:11]=[C:13]([NH2:15])[CH:16]=[CH:14][CH:12]1)[C:7](=[O:9])[OH:10]).[Si:4]([CH3:6])[CH3:8].[ClH:3]";
		//String reactionFile = "[cH2:1]1([Br:21])[cH2:2][cH2:3][cH:4]([cH2:5][cH2:6]1)[cH:7]2[cH2:8][cH2:9][cH2:10]([Br:22])[cH2:11][cH2:12]2.[NH2:13][cH:14]3[cH2:15][cH2:16][cH2:17][cH:18]([cH2:19]3)[CH3:20]>>[cH2:1]1([NH2:13][cH:14]3[cH2:15][cH2:16][cH2:17][cH:18]([cH2:19]3)[CH3:20])[cH2:2][cH2:3][cH:4]([cH2:5][cH2:6]1)[cH:7]2[cH2:8][cH2:9][cH2:10]([NH2:13][cH:14]3[cH2:15][cH2:16][cH2:17][cH:18]([cH2:19]3)[CH3:20])[cH2:11][cH2:12]2";
		//String reactionFile = "/Users/delanneevc/Downloads/valid.txt";
		String reactionFile = "[CH3:1][C:2](=[O:3])[CH3:4].[CH3:10][Si:4][O:5][S:6](=[O:7])(=[O:8])[CH3:9]>>[CH2:1]=[C:2]([O:3][Si:4][CH3:10])[CH3:4]";
		//2 reactants -> 3 products but one missing (generated a wrong PM because of RC detection) CORRECTED
		//String reactionFile = "[CH3:1][Si:4]([Cl:3])([CH3:6])[CH3:8].[CH2:2]([C:5]=1[CH:11]=[C:13]([NH2:15])[CH:16]=[CH:14][CH:12]1)[C:7](=[O:9])[OH:10]>>[CH2:1]([CH:2]([C:5]=1[CH:11]=[C:13]([NH2:15])[CH:16]=[CH:14][CH:12]1)[C:7](=[O:9])[OH:10]).[Cl:3]";
		//test substitution case
		//String outputDirectory = "/datahdd/SAVI/Alexey/testMapping/";
		String outputDirectory = "/Users/delanneevc/Documents/testReactionCode/";
		String errorsFileName;
		boolean makePseudoSmiles = false;
		boolean makeImage = false;
		boolean onlyPseudoMol = false;
		boolean writeSDF = false;
		String prefix = "PR";
		String outputFormat = "csv";


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

		FileWriter outputfile;
		CSVWriter csvWriter = null;
		
		ObjectMapper mapper = new ObjectMapper();
		List<Object> results = new ArrayList<Object>();
		List<String> smiles = new ArrayList<String>();
		
		if (outputFormat.equals("csv")) {
			outputfile = new FileWriter(outputDirectory + prefix + "_reactionsCode.csv"); 
			System.out.println(outputDirectory + prefix + "_reactionsCode.csv");
			csvWriter = new CSVWriter(outputfile, '\t', CSVWriter.DEFAULT_QUOTE_CHARACTER, CSVWriter.DEFAULT_ESCAPE_CHARACTER,
					CSVWriter.DEFAULT_LINE_END); 
		}
		else if (outputFormat.equals("json")) {
			outputfile = new FileWriter(outputDirectory + prefix + "_reactionsCode.json"); 
			
		}

		List<String> headers = new ArrayList<String>();
		
		//parse reactions
		ChemicalFormatParser parser = new ChemicalFormatParser();
		IReactionSet reactions = DefaultChemObjectBuilder.getInstance().newInstance(IReactionSet.class);
		String format = parser.formatDetector(reactionFile);

		if (format.equals("RDF_FILE")) {
			reactions = parser.parseRDF(reactionFile);
		}
		else if (format.equals("RXN")) {
			reactions = parser.parseRXN(reactionFile);
		}
		else if (format.equals("SMIRKS_FILE")) {
			reactions = parser.parseReactionSMILES(reactionFile, 0, true, false);
		}
		else if (format.equals("SMIRKS_TEXT")) {
			reactions = parser.parseReactionSMILES(reactionFile, false);
		}
		else {
			System.err.println("Wrong input file: file has to be formated as an RXN or an RDF file "+format);
			return;
		}

		for (int i = 0; i < reactions.getReactionCount(); i++) {
			if (i == 500) {
				break;
			}
			IReaction reaction = reactions.getReaction(i);
//			System.out.println("reaction " + i );
//			System.out.println(reaction.getProperties());
			Map<Object, Object> properties = reaction.getProperties();
			IAtomContainerSet reactants = reaction.getReactants();
			IAtomContainerSet products = reaction.getProducts();

			for (IAtomContainer ac : reactants.atomContainers()) {
//				tools.perceiveAromaticityAndflagAtomsAndBonds(ac);
				tools.attributeIDtoAtomsAndBonds(ac);
			}
			for (IAtomContainer ac : products.atomContainers()) {
//				tools.perceiveAromaticityAndflagAtomsAndBonds(ac);
				tools.attributeIDtoAtomsAndBonds(ac);
			}
			//Annotate the reaction with right constants
			PseudoMolecule pm = new PseudoMolecule();
			IAtomContainer pseudoMolecule = null;

			//try {
				pm.reactionAnnotator(reaction);
				pseudoMolecule = pm.makePseudoMolecule(reactants, products);
				pseudoMolecule.setProperty(CDKConstants.TITLE, "PSEUDO-MOLECULE");
				
			//}
			//catch (Exception e) {
				//dataFile(reaction.getID() + "\t Something get wrong during the PseudoMoleclule generation\n", errorsFileName);
				//continue;
			//}

			if (writeSDF == true) {
				try {
					pm.writePseudoMoleculeSDFile(prefix+"_"+i, new File(outputDirectory, "SDF_PseudoMolecule").toString(), 
							pseudoMolecule);
				}
				catch (Exception e) {
					//dataFile(reaction.getID() + "\t The SDF could not be generated\n", errorsFileName);
				}
			}
			if (makeImage == true) {
				try {
					pm.writePseudoMoleculeImage(new File(outputDirectory, "IMG_PseudoMolecule").toString(), prefix+"_"+i);
				}
				catch (Exception e) {
					System.out.println("err");
//					dataFile(reaction.getID() + "\t The IMG could not be generated\n", errorsFileName);
				}
			}

			String pseudoSmiles = null;
			try {
				pseudoSmiles = pm.GetPseudoMoleculeSmiles(pseudoMolecule);
			}
			catch (Exception e) {
				//						dataFile(reaction.getID() + "\t The pseudoSmiles could not be generated\n", errorsFileName);
			}

			if (makePseudoSmiles) {
				if (pseudoSmiles != null) {
					smiles.add(pseudoSmiles);
				}
			}
			if (onlyPseudoMol == false) {
				//get atom duplication number
				Map<String, Integer> numberOfRepetitions = pm.atomRepetition();
				Set<IAtom> reactionCenterAtom = pm.getReactioncenter();
				
				if (reactionCenterAtom.size() == 0) {
					System.out.println(reaction.getID() + " No atom has been identified as a reaction centre. The mapping is probably incorrect\n");
					System.out.println(reaction.getProperties());
					break;
					//							dataFile(reaction.getID() + "\t No atom has been identified as a reaction centre. The mapping is probably incorrect\n", errorsFileName);
					//continue;
				}
//				try {
				//set parameters for generation of reactionCode
				EncodeReactionCode encoder = new EncodeReactionCode();
				encoder.setBondType(true);
				encoder.setCharge(true);
				encoder.setHybridization(false);
				encoder.setRepetition(true);
				encoder.setStereochemistry(true);
					//make reactionCode
					Map<String,String> reactionCodeMap = encoder.generateReactionCode(reactionCenterAtom, reactants, 
							pseudoMolecule, numberOfRepetitions);
					//String reactionCode = ReactionCode.reactionCodeMapToStringMultiLines(reactionCodeMap);
					String reactionCode = encoder.reactionCodeMapToStringOneLine(reactionCodeMap);
					System.out.println(reactionCode);
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
							if (i == 0) {
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
//				}
//				catch (Exception e) {
//					dataFile(reaction.getID() + "\t Something get wrong during the ReactionCode generation\n", errorsFileName);
//					continue;
//				}


			}
		}

		if (outputFormat.equals("json")) {
			//write JSON
			try (FileWriter file = new FileWriter(outputDirectory + prefix + "_reactionsCode.json")) {
				mapper.writerWithDefaultPrettyPrinter().writeValue(file, results);
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
		else if (outputFormat.equals("csv")) {
			csvWriter.close();
		}

		if (makePseudoSmiles) {
			FileWriter writer = new FileWriter(outputDirectory + prefix + "_pseudoSmiles.rxt"); 
			for(String str : smiles) {
				writer.write(str+"\n");
			}
			writer.close();
		}
	}
}

