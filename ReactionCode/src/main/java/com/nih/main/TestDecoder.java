package com.nih.main;


import java.awt.Color;
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
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.graph.ConnectivityChecker;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IReaction;
import org.openscience.cdk.interfaces.IReactionSet;

import com.nih.codes.DecodeReactionCode;
import com.nih.codes.EncodeReactionCode;
import com.nih.parser.ChemicalFormatParser;
import com.nih.reaction.PseudoMolecule;
import com.nih.tools.tools;
import com.nih.transformer.Transformer;
import com.nih.writer.MDLV2000RDFWriter;
import com.nih.writer.MDLV2000RXNWriter;
import com.opencsv.CSVWriter;

import uk.ac.ebi.beam.Atom;
import uk.ac.ebi.beam.Bond;
import uk.ac.ebi.beam.Edge;
import uk.ac.ebi.beam.Graph;

public class TestDecoder {

	public static void main(String[] args) throws IOException, CDKException, CloneNotSupportedException {
		// TODO Auto-generated method stub

//		String reactionFile = "/datahdd/SAVI/Alexey/testMapping/testMapping.rdf";
//		String reactionFile = "/home/delanneev/Downloads/reaxys_10038803_exampleForReactionCode.rdf";
//		String reactionFile = "/home/delanneev/Downloads/reaxys_10038803_exampleForReactionCode_WithIsotopeAndCharge.rdf";
//		String reactionFile = "/home/delanneev/Downloads/test1.rxn";
		//String reactionFile = "﻿[CH3:26][N:27]1[CH2:28][CH2:29][CH2:30][C:31]1=[O:32].[Cl:1][c:2]1[n:3](-[c:12]2[cH:13][cH:14][c:15]([Cl:18])[cH:16][cH:17]2)[n:4][c:5]2[cH:6][c:7]([F:11])[cH:8][cH:9][c:10]12.[NH2:19][CH:20]1[CH2:21][CH2:22][CH2:23][CH2:24][CH2:25]1>>[c:2]1([NH:19][CH:20]2[CH2:21][CH2:22][CH2:23][CH2:24][CH2:25]2)[n:3](-[c:12]2[cH:13][cH:14][c:15]([Cl:18])[cH:16][cH:17]2)[n:4][c:5]2[cH:6][c:7]([F:11])[cH:8][cH:9][c:10]12";
//		String reactionFile = "[CH2:42]([N:43]=[C:44]=[N:45][CH2:46][CH2:47][CH2:48][N:49]([CH3:50])[CH3:51])[CH3:52].[CH3:32][N:33]([CH2:34][CH:35]=[CH:36][C:37](=[O:38])[OH:39])[CH3:40].[CH3:64][N:65]([CH3:66])[CH:67]=[O:68].[CH3:70][CH2:71][N:72]([CH2:73][CH3:74])[CH2:75][CH3:76].[ClH:1].[ClH:2].[ClH:31].[ClH:41].[NH2:3][CH2:4][c:5]1[cH:6][c:7]2[n:8][cH:9][n:10][c:11]([NH:14][c:15]3[cH:16][c:17]([Cl:30])[c:18]([O:21][CH2:22][c:23]4[cH:24][c:25]([F:29])[cH:26][cH:27][cH:28]4)[cH:19][cH:20]3)[c:12]2[nH:13]1.[OH2:53].[OH2:69].[OH:54][n:55]1[c:56]2[cH:57][cH:58][cH:59][cH:60][c:61]2[n:62][n:63]1>>[NH:3]([CH2:4][c:5]1[cH:6][c:7]2[n:8][cH:9][n:10][c:11]([NH:14][c:15]3[cH:16][c:17]([Cl:30])[c:18]([O:21][CH2:22][c:23]4[cH:24][c:25]([F:29])[cH:26][cH:27][cH:28]4)[cH:19][cH:20]3)[c:12]2[nH:13]1)[C:37]([CH:36]=[CH:35][CH2:34][N:33]([CH3:32])[CH3:40])=[O:38]";
//		String reactionFile = "[C:1]([c:2]1[cH:3][cH:4][cH:5][cH:6][cH:7]1)(=[O:8])[Cl:9].[CH3:10][CH2:11][CH2:12][NH2:13].[Cl:14][CH2:15][Cl:16]>>[C:1]([c:2]1[cH:3][cH:4][cH:5][cH:6][cH:7]1)(=[O:8])[NH:13][CH2:12][CH2:11][CH3:10]";
//		String reactionFile = "[CH2:23]1[O:24][CH2:25][CH2:26][CH2:27]1.[F:1][c:2]1[c:3]([N+:10](=[O:11])[O-:12])[cH:4][c:5]([F:9])[c:6]([F:8])[cH:7]1.[H-:22].[NH2:13][c:14]1[s:15][cH:16][cH:17][c:18]1[C:19]#[N:20].[Na+:21]>>[c:2]1([NH:13][c:14]2[s:15][cH:16][cH:17][c:18]2[C:19]#[N:20])[c:3]([N+:10](=[O:11])[O-:12])[cH:4][c:5]([F:9])[c:6]([F:8])[cH:7]1";
		String reactionFile = "[Br:25][c:26]1[cH:27][cH:28][c:29]2[cH:30][cH:31][cH:32][c:33]3[c:34]4[cH:35][cH:36][cH:37][c:38]5[cH:39][cH:40][cH:41][c:42]([c:43]1[c:44]23)[c:45]45.[cH:1]1[cH:2][c:3]2[cH:4][c:5]3[cH:6][cH:7][c:8]([cH:9][c:10]4[cH:11][cH:12][c:13]([cH:14][c:15]5[cH:16][cH:17][c:18]([cH:19][c:20]1[n:21]2)[nH:22]5)[n:23]4)[nH:24]3>>[cH:1]1[cH:2][c:3]2[cH:4][c:5]3[cH:6][cH:7][c:8]([cH:9][c:10]4[cH:11][cH:12][c:13]([cH:14][c:15]5[cH:16][cH:17][c:18]([cH:19][c:20]1[nH:21]2)[n:22]5)[nH:23]4)[n:24]3.[cH:26]1[cH:27][cH:28][c:29]2[cH:30][cH:31][cH:32][c:33]3[c:34]4[cH:35][cH:36][cH:37][c:38]5[cH:39][cH:40][cH:41][c:42]([c:43]1[c:44]23)[c:45]45";
		String outputDirectory = "/datahdd/SAVI/Alexey/testMapping/Decoder/";
		String errorsFileName;
		boolean makeImage = true;
		boolean writeSMIRKS = true;
		boolean writeSMIRKSWMapping = true;
		boolean writeRXN = true;
		boolean writeRDF = true;
		String prefix = "PR";
		
		reactionFile = reactionFile.replaceAll("[^\\x00-\\x7F]", "");
		
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
			directory = new File(outputDirectory, "IMG_OriginalReaction");
			if (!directory.exists()){
				directory.mkdir();
			}
		}

		List<String> reactionCodes = new ArrayList<String>();

		System.out.println("Encode Reaction");
		//parse reactions
		ChemicalFormatParser parser = new ChemicalFormatParser();
		IReactionSet reactions = DefaultChemObjectBuilder.getInstance().newInstance(IReactionSet.class);
		String format = parser.formatDetector(reactionFile);
		
		if (format.equals("RDF_FILE")) {
			reactions = parser.parseRDF(reactionFile);
		}
		else if (format.equals("RXN_FILE")) {
			reactions = parser.parseRXN(reactionFile);
		}
		else if (format.equals("SMIRKS_TEXT")) {
			reactions = parser.parseReactionSMILES(reactionFile, true);
		}
		else {
			System.err.println("Wrong input file: file has to be formated as an RXN or an RDF file");
			return;
		}

		for (int i = 0; i < reactions.getReactionCount(); i++) {
			if (i == 500) {
				break;
			}
			
			IReaction reaction = reactions.getReaction(i);
			IAtomContainerSet reactants = reaction.getReactants();
			IAtomContainerSet products = reaction.getProducts();
//			
//			for (IAtomContainer ac : products.atomContainers()) {
//				for (IAtom a : ac.atoms()) {
//					if ((int) a.getProperty(CDKConstants.ATOM_ATOM_MAPPING) == 4) {
//						ac.removeBond(ac.getConnectedBondsList(a).get(0));
//						ac.removeAtom(a);
//						break;
//					}
//				}
//			}
			
			if (makeImage) {
				new org.openscience.cdk.depict.DepictionGenerator()
				.withOuterGlowHighlight().withAtomColors().withAtomMapNumbers()
				.depict(reaction).writeTo(outputDirectory + "/IMG_OriginalReaction/reaction_" + i + ".pdf");
			}

			for (IAtomContainer ac : reactants.atomContainers()) {
				//tools.perceiveAromaticityAndflagAtomsAndBonds(ac);
				tools.attributeIDtoAtomsAndBonds(ac);
			}
			for (IAtomContainer ac : products.atomContainers()) {
				//tools.perceiveAromaticityAndflagAtomsAndBonds(ac);
				tools.attributeIDtoAtomsAndBonds(ac);
			}
			//Annotate the reaction with right constants
			PseudoMolecule pm = new PseudoMolecule();
			IAtomContainer pseudoMolecule = null;

			try {
				pm.reactionAnnotator(reaction);
				pseudoMolecule = pm.makePseudoMolecule(reactants, products);
				pseudoMolecule.setProperty(CDKConstants.TITLE, "PSEUDO-MOLECULE");
				//pm.writePseudoMoleculeImage(new File(outputDirectory).toString(), prefix+"_"+"000");
			}
			catch (Exception e) {
				//dataFile(reaction.getID() + "\t Something get wrong during the PseudoMoleclule generation\n", errorsFileName);
				continue;
			}


			//get atom duplication number
			Map<String, Integer> numberOfRepetitions = pm.atomRepetition();
			Set<IAtom> reactionCenterAtom = pm.getReactioncenter();

			if (reactionCenterAtom.size() == 0) {
				continue;
			}
			//						try {
			//make reactionCode
			EncodeReactionCode encoder = new EncodeReactionCode();
			Map<String,String> reactionCodeMap = encoder.makeReactionCode(reactionCenterAtom, reactants, 
					pseudoMolecule, numberOfRepetitions);
			//String reactionCode = ReactionCode.reactionCodeMapToStringMultiLines(reactionCodeMap);
			String reactionCode = encoder.reactionCodeMapToStringOneLine(reactionCodeMap);
//			reactionCodes.add(reactionCode);
			System.out.println(reactionCode);
		}

		//parse reactions
		System.out.println("Decode ReactionCode");

		//		IReactionSet reactions = parser.parseReactionCode(reactionFile);

		//String reactionCode = "0:906()[1]708(11GG)[1]508(01GG)[1]506(21GI10GH)[1]|1:006(11GJ)[1]006(11GG)[1]|2:006(11GL11GK)[1]|3:006(11GM)[1]006(11GM)[1]|4:006(11GO11GN)[1]|"; 
		//String reactionCode = "0:907()[1]708()[1]706(10GH01GG)[1]|1:008(22GI)[1]006(11GI)[1]|2:006(99GK)[1]006(99GK)[1]|3:007(11GM)[1]006(99GL)[1]006(99GM)[1]/c0-II;|4:008(22GN)[1]008(11GN)[1]006(99GP99GO)[1]/c2-HH;|5:008(11GS)[1]|6:006(11GT)[1]|A:006(11GG)[1]006(11GG)[1]006(11GG)[1]|B:008(2200)[1]|";
//		String reactionCode = "0:906()[1]723()[1]706(10GH01GG)[1]|1:007(92GI)[1]006(91GI)[1]006(22GG)[1]";
		//String reactionCode = "0:906()[1]723()[1]706(10GH01GG)[1]|1:007(92GI)[1]006(91GI)[1]006(22GG)[1]|2:006(92GK)[1]006(91GJ)[1]006(11GL)[1]/s5-08;|3:008(22GO)[1]008(11GO)[1]006(92GN91GM)[1]006(11GN)[1]|4:008(11GS)[1]006(11GQ)[1]|5:006(11GU)[1]|";
		//String reactionCode = "0:906()[1]906(01GG)[1]|1:006(11GH)[1]006(11GH)[1]|2:011(1100)[1]008(11GJ)[1]008(22GJ)[1]006(99GI)[1]006(99GI)[1]|";
		//String reactionCode = "0:906()[1]906(01GG)[1]|1:006(11GH)[1]006(11GH)[1]|2:011(1100)[1]008(11GJ)[1]008(22GJ)[1]006(99GI)[1]006(99GI)[1]|3:006(99GN)[1]006(99GO)[1]|4:007(11GQ)[1]006(99GQ99GP)[1]|A:00E(11GG)[1]|B:006(1100)[1]006(1100)[1]|";
		//String reactionCode = "0:906()[1]711()[1]70E(10GH)[1]706(10GI01GG)[1]|1:006(11GG)[1]006(11GI)[1]006(11GG)[1]006(11GI)[1]|2:008(22GM)[1]008(11GM)[1]006(99GK)[1]006(99GK)[1]|3:006(99GR)[1]006(99GQ)[1]|4:007(11GS)[1]006(99GT99GS)[1]|";
		//String reactionCode = "0:906()[1]708()[1]708()[1]707()[1]707(99GJ)[1]707(01GG)[1]706(90GJ)[1]706(90GK)[1]706()[1]706()[1]706(90GP90GO)[1]706(11GM10GL)[1]706()[1]706(02GP01GS01GO)[1]706(02GK01GS)[1]506(11GG10GS10GH02GJ)[1]506(11GV10GU10GT10GI02GS)[1]|1:009(11GP)[1]006(92GO)[1]006(91GP)[1]006(91GL)[1]006(91GL)[1]|2:011(1101)[1]007(92HG)[1]006(92GZ91GY)[1]006(92HH91HJ)[1]006(11GZ)[1]006(11HG)[1]|3:009(11HM)[1]009(11HM)[1]009(11HM)[1]|A:006(99GM)[1]006(99GN9900)[1]005(11GQ11GI11GH)[1]|";
		//String reactionCode = "0:907()[1]708()[1]706(10GH01GG)[1]|1:008(22GI)[1]006(11GG)[1]|2:006(99GK)[1]006(99GK)[1]|3:006(99GL)[1]006(99GM)[1]|4:006(99GO99GN)[1]|5:006(11GP)[1]|6:006(99GQ)[1]006(99GQ)[1]|7:009(11GR)[1]006(99GR)[1]006(99GS)[1]|8:007(91GV)[1]006(99GV99GU)[1]|9:006(91GW)[1]006(91GX)[1]006(11GW)[1]|10:008(22GZ)[1]006(92GY91GZ)[1]006(11HG)[1]|11:006(11HI)[1]|12:008(22HK)[1]008(11HK)[1]|";
		//String reactionCode = "0:907()[1]709()[1]706(10GH01GG)[1]|1:006(99GI)[1]006(99GI)[1]006(11GG)[1]|2:010(99GL)[1]007(11GJ)[1]006(99GL)[1]006(99GJ)[1]006(99GK)[1]/c2-II;|3:009(11GQ)[1]008(22GN)[1]008(11GN)[1]006(99GM)[1]006(99GU99GO)[1]006(99GQ99GP)[1]006(11GO)[1]/c4-HH;|4:009(11GW)[1]007(33GX)[1]|";
		//String reactionCode = "0:907()[1]711()[1]706(10GH01GG)[1]|1:007(99GI)[1]006(99GI)[1]006(11GG)[1]|2:007(99GJ)[1]006(99GM99GK)[1]006(11GJ)[1]006(99GK)[1]006(11GL)[1]006(11GL)[1]|3:006(99GN)[1]006(99GO)[1]006(99GO)[1]006(99GP)[1]006(11GQ)[1]006(11GR)[1]|4:006(99GV99GS)[1]006(99GT)[1]006(99GU)[1]006(11GX11GW)[1]|5:009(11GY)[1]006(99HG99GZ)[1]|6:011(11HJ)[1]|";
		//String reactionCode = "0:723()[1]706(10GG)[1]|1:006(99GH)[1]006(99GH)[1]|2:006(99GI)[1]006(99GI)[1]006(99GJ)[1]|3:006(99GM99GL)[1]006(99GK)[1]006(99GK)[1]006(99GL)[1]|4:006(99GQ99GO)[1]006(99GN)[1]006(99GO)[1]006(99GP)[1]006(99GQ)[1]|5:006(99GV99GS)[1]006(99GU99GT)[1]006(99GR)[1]006(99GT)[1]|6:006(99GZ99GY)[1]|";
		String reactionCode = "0:70E()[1]708(10GG)[1]508(01GG)[1]506(21GI)[1]506(12GJ)[1]|1:006(11GJ)[1]006(11GG)[1]006(11GG)[1]006(11GG)[1]/i06II|2:008(11GM)[1]008(11GM)[1]006(22GL)[1]006(11GL)[1]006(11GM)[1]/c00HH/i00JJ02HH|3:006(11GS)[1]006(11GU11GR)[1]/s00210212|4:008(11GU)[1]006(11GV)[1]/s01640364|5:006(11GW)[1]006(22GX)[1]006(11GX)[1]|A:010(11GH)[1]|B:008(2200)[1]008(2200)[1]006(1100)[1]|C:009(1103)[1]009(1103)[1]009(1103)[1]|";
		reactionCodes.add(reactionCode);
		DecodeReactionCode decoder = new DecodeReactionCode();
		IReactionSet reactions2 = decoder.makeReactions(reactionCodes);	
		
//		Transformer transformer = new Transformer();
//		IReaction ptrn = decoder.decode("0:907()[1]711()[1]706(10GH01GG)[1]|1:007(99GI)[1]006(99GI)[1]006(11GG)[1]");
//		IReaction trans = transformer.applyTranform("[NH2:1][CH:6]1[CH2:11][CH2:17][CH2:22][CH2:18][CH2:12]1.[Cl:2][cH:3]1[n:4]([nH:7][cH:8]2[cH:5]1[cH2:10][cH2:16][cH:19]([cH2:13]2)[F:23])-[cH:9]3[cH2:14][cH2:20][cH:24]([cH2:21][cH2:15]3)[Cl:25]", ptrn);
//		System.out.println(tools.makeSmiles(trans, false));

		if (writeRDF) {
			File file = new File(outputDirectory + "decodedReactions.rdf");
			try (MDLV2000RDFWriter writer = new MDLV2000RDFWriter(new FileWriter(file))) {
				writer.setWriteAromaticBondTypes(true);
				writer.write(reactions2);
				//writer.setRdFieldsProperties(map);
				writer.close();
			} catch (CDKException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		if (writeRXN || makeImage || writeSMIRKS || writeSMIRKSWMapping) {
			FileWriter smirksWriter = null;
			FileWriter smirksWriterWithMapping = null;
			if (writeSMIRKS)
				smirksWriter = new FileWriter(outputDirectory + prefix + "_reactionSMILES.txt"); 
			if (writeSMIRKSWMapping)
				smirksWriterWithMapping = new FileWriter(outputDirectory + prefix + "_reactionSMILESAndMapping.txt"); 
			for (int i = 0; i < reactions2.getReactionCount(); i++) {
				IReaction reaction = reactions2.getReaction(i);
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
					new org.openscience.cdk.depict.DepictionGenerator()
					.withHighlight(decoder.getBondsFormed(), Color.BLUE)
					.withHighlight(decoder.getBondsCleaved(), Color.RED)
					.withHighlight(decoder.getBondsOrder(), Color.GREEN)
					.withHighlight(decoder.getReactionCenter(), Color.LIGHT_GRAY)
					.withOuterGlowHighlight().withAtomColors().withAtomMapNumbers()
					.depict(reaction).writeTo(outputDirectory + "/IMG_ReactionDecoded/reaction_" + i + ".pdf");
				}
				if (writeSMIRKSWMapping) {
					smirksWriterWithMapping.write(tools.makeSmiles(reaction, true) + "\n");
					System.out.println(tools.makeSmiles(reaction, true));
				}
				if (writeSMIRKS) {
					smirksWriter.write(tools.makeSmiles(reaction, false) + "\n");
				}
			}
			if (writeSMIRKS)
				smirksWriter.close(); 
			if (writeSMIRKSWMapping)
				smirksWriterWithMapping.close(); 
		}
	}
	
	
}
