package com.nih.main;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.List;
import java.util.Map;

import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.aromaticity.Aromaticity;
import org.openscience.cdk.aromaticity.ElectronDonation;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.fingerprint.CircularFingerprinter;
import org.openscience.cdk.fingerprint.GraphOnlyFingerprinter;
import org.openscience.cdk.fingerprint.IBitFingerprint;
import org.openscience.cdk.fingerprint.SignatureFingerprinter;
import org.openscience.cdk.graph.CycleFinder;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IReaction;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.smiles2.SmilesParser2;

import com.nih.parser.ChemicalFormatParser;
import com.nih.reaction.PseudoMolecule;
import com.nih.reaction.additionalConstants;
import com.nih.tools.tools;

import uk.ac.ebi.beam.Graph;

public class Tests {
	
	static List<Integer> err = new ArrayList<Integer>();
	
	public static void main(String[] args) throws IOException, CDKException, CloneNotSupportedException {
		String reactionFile = "/Users/delanneevc/Downloads/USPTO_SMIRKS.txt";
		
		//test1(reactionFile);
		test2(reactionFile);
		//test3(reactionFile);
		
	}
	
	
	private static void test1(String reactionFile) throws CDKException, IOException, CloneNotSupportedException {
		SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
		
		BufferedReader br = new BufferedReader(new FileReader(reactionFile)); 
		String line;
		int error = 0;
		int success = 0;
		int cpt = 0;
		
		while ((line = br.readLine()) != null) {
			IReaction reaction1 = sp.parseReactionSmiles(line);
			API api = new API();
			String reactionCode = api.encode((IReaction) reaction1.clone(), true, true, false, true, true, false);
			IReaction reaction2 = api.decode(reactionCode, false, false);
			
			//remove aromaticity
			removeAromaticity(reaction1);
			removeAromaticity(reaction2);
			
			removeAgents(reaction1, line.split(" ")[1]);
			
			List<BitSet> l1 = new ArrayList<BitSet>();
			List<BitSet> l2 = new ArrayList<BitSet>();
			
			for (IAtomContainer ac : reaction1.getReactants().atomContainers()) {
				CircularFingerprinter fp = new CircularFingerprinter();
				fp.calculate(ac);
				l1.add(fp.getBitFingerprint(ac).asBitSet());
			}
			for (IAtomContainer ac : reaction1.getProducts().atomContainers()) {
				CircularFingerprinter fp = new CircularFingerprinter();
				fp.calculate(ac);
				l1.add(fp.getBitFingerprint(ac).asBitSet());
			}
			for (IAtomContainer ac : reaction2.getReactants().atomContainers()) {
				CircularFingerprinter fp = new CircularFingerprinter();
				fp.calculate(ac);
				l2.add(fp.getBitFingerprint(ac).asBitSet());
			}
			for (IAtomContainer ac : reaction2.getProducts().atomContainers()) {
				CircularFingerprinter fp = new CircularFingerprinter();
				fp.calculate(ac);
				l2.add(fp.getBitFingerprint(ac).asBitSet());
			}
			l1.removeAll(l2);
			if (!l1.isEmpty()) {
				err.add(cpt);
				System.out.println("different reactions "+cpt);
				System.out.println(reactionCode);
				System.out.println(line);
				System.out.println(tools.makeSmiles(reaction1, true));
				System.out.println(" "+tools.makeSmiles(reaction2, true));
				error++;
			}
			else {
				success++;
			}
			cpt++;
		}
		System.out.println("valids: " + success);
		System.out.println("errors: " + error);
	}
	
	private static void test2(String reactionFile) throws CDKException, IOException, CloneNotSupportedException {
		SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
		
		BufferedReader br = new BufferedReader(new FileReader(reactionFile)); 
		String line;
		int error = 0;
		int kerror = 0;
		int success = 0;
		int cpt = 0;
		
		while ((line = br.readLine()) != null) {
			IReaction reaction1 = sp.parseReactionSmiles(line);
			API api = new API();
			String reactionCode = api.encode((IReaction) reaction1.clone(), true, true, false, true, true, false);
			IReaction reaction2 = api.decode(reactionCode, true, false);
			
			removeAgents(reaction1, line.split(" ")[1]);
			
			List<BitSet> l1 = new ArrayList<BitSet>();
			List<BitSet> l2 = new ArrayList<BitSet>();
			
			for (IAtomContainer ac : reaction1.getReactants().atomContainers()) {
				CircularFingerprinter fp = new CircularFingerprinter();
				fp.calculate(ac);
				l1.add(fp.getBitFingerprint(ac).asBitSet());
			}
			for (IAtomContainer ac : reaction1.getProducts().atomContainers()) {
				CircularFingerprinter fp = new CircularFingerprinter();
				fp.calculate(ac);
				l1.add(fp.getBitFingerprint(ac).asBitSet());
			}
			for (IAtomContainer ac : reaction2.getReactants().atomContainers()) {
				CircularFingerprinter fp = new CircularFingerprinter();
				fp.calculate(ac);
				l2.add(fp.getBitFingerprint(ac).asBitSet());
			}
			for (IAtomContainer ac : reaction2.getProducts().atomContainers()) {
				CircularFingerprinter fp = new CircularFingerprinter();
				fp.calculate(ac);
				l2.add(fp.getBitFingerprint(ac).asBitSet());
			}
			l1.removeAll(l2);
			if (!l1.isEmpty()) {
				System.out.println("different reactions "+cpt);
				System.out.println(reactionCode);
				System.out.println(line);
				System.out.println(tools.makeSmiles(reaction1, true));
				String decodedSMIRKS = tools.makeSmiles(reaction2, true);
				if (!isKekulizationValid(decodedSMIRKS))
					kerror++;
				System.out.println(" "+decodedSMIRKS);
				error++;
			}
			else {
				//System.out.println(" "+tools.makeSmiles(reaction2, true));
				success++;
			}
			cpt++;
		}
		int tauto = error - kerror;
		System.out.println("valids: " + success);
		System.out.println("errors: " + error);
		System.out.println("errors related unsuccesful kekulization: " + kerror);
		System.out.println("errors related to different kekule: " + tauto);
	}
	
	private static void test3(String reactionFile) throws CDKException, IOException, CloneNotSupportedException {
		SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
		
		BufferedReader br = new BufferedReader(new FileReader(reactionFile)); 
		String line;
		int error = 0;
		int success = 0;
		int cpt = 0;
		
		while ((line = br.readLine()) != null) {
			IReaction reaction1 = sp.parseReactionSmiles(line);
			API api = new API();
			String reactionCode = api.encodeAndRebalance((IReaction) reaction1.clone(), true, true, false, true, true, false);
			IReaction reaction2 = api.decode(reactionCode, false, true);
			//remove aromaticity
			removeAromaticity(reaction1);
			removeAromaticity(reaction2);
			
			removeAgents(reaction1, line.split(" ")[1]);
			
			List<BitSet> l1 = new ArrayList<BitSet>();
			List<BitSet> l2 = new ArrayList<BitSet>();
			
			for (IAtomContainer ac : reaction1.getReactants().atomContainers()) {
				CircularFingerprinter fp = new CircularFingerprinter();
				fp.calculate(ac);
				l1.add(fp.getBitFingerprint(ac).asBitSet());
			}
			for (IAtomContainer ac : reaction1.getProducts().atomContainers()) {
				CircularFingerprinter fp = new CircularFingerprinter();
				fp.calculate(ac);
				l1.add(fp.getBitFingerprint(ac).asBitSet());
			}
			for (IAtomContainer ac : reaction2.getReactants().atomContainers()) {
				CircularFingerprinter fp = new CircularFingerprinter();
				fp.calculate(ac);
				l2.add(fp.getBitFingerprint(ac).asBitSet());
			}
			for (IAtomContainer ac : reaction2.getProducts().atomContainers()) {
				CircularFingerprinter fp = new CircularFingerprinter();
				fp.calculate(ac);
				l2.add(fp.getBitFingerprint(ac).asBitSet());
			}
			l1.removeAll(l2);
			if (!l1.isEmpty()) {
				if (!err.contains(cpt)) {
					System.out.println("different reactions "+cpt);
					System.out.println(reactionCode);
					System.out.println(line);
					System.out.println(tools.makeSmiles(reaction1, true));
					System.out.println(" "+tools.makeSmiles(reaction2, true));
				}
				error++;
			}
			else {
				success++;
			}
			cpt++;
		}
		System.out.println("valids: " + success);
		System.out.println("errors: " + error);
	}
	
	private static void removeAromaticity(IReaction reaction) {
		for (IAtomContainer ac : reaction.getReactants().atomContainers()) {
			removeAromaticity(ac);
		}
		for (IAtomContainer ac : reaction.getProducts().atomContainers()) {
			removeAromaticity(ac);
		}
	}
	
	private static void removeAromaticity(IAtomContainer ac) {
		for (IBond bond : ac.bonds()) {
			if (bond.isAromatic()) {
				bond.setOrder(IBond.Order.SINGLE);
				bond.setIsAromatic(false);
				bond.getBegin().setImplicitHydrogenCount(0);
				bond.getEnd().setImplicitHydrogenCount(0);
				bond.getBegin().setIsAromatic(false);
				bond.getEnd().setIsAromatic(false);
			}
		}
	}
	
	private static void removeAgents(IReaction reaction, String rc) {
		List<Integer> aam = new ArrayList<Integer>();
		for (String s1 : rc.split(";")) {
			for (String aam2 : s1.split("-")) {
				aam.add(Integer.parseInt(aam2));
			}
		}
		List<Integer> aamInProducts = new ArrayList<Integer>();
		List<IAtomContainer> toRemove = new ArrayList<IAtomContainer>();
		List<IAtomContainer> toRemove2 = new ArrayList<IAtomContainer>();
		for (IAtomContainer ac : reaction.getProducts().atomContainers()) {
			boolean isPartOfAReactionCenter = false;
			for (IAtom a : ac.atoms()) {
				int map = a.getProperty(CDKConstants.ATOM_ATOM_MAPPING);
				aamInProducts.add(map);
				if (aam.contains(map))
					isPartOfAReactionCenter = true;
			}
			if (!isPartOfAReactionCenter)
				toRemove2.add(ac);
		}
		for (IAtomContainer ac : reaction.getReactants().atomContainers()) {
			boolean isAgent = true;
			boolean isPartOfAReactionCenter = false;
			for (IAtom a : ac.atoms()) {
				int map = a.getProperty(CDKConstants.ATOM_ATOM_MAPPING);
				if (aamInProducts.contains(map))
					isAgent = false;
				if (aam.contains(map))
					isPartOfAReactionCenter = true;
			}
			if (isAgent || !isPartOfAReactionCenter)
				toRemove.add(ac);
		}
		for (IAtomContainer ac : toRemove) {
			reaction.getReactants().removeAtomContainer(ac);
		}
		for (IAtomContainer ac : toRemove2) {
			reaction.getProducts().removeAtomContainer(ac);
		}
	}
	
	private static boolean isKekulizationValid(String smirks) {
		SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
		IReaction reaction;
		// will return false if the kekule is false as CDK will fail to parse the SMIRKS
		try {
			reaction = sp.parseReactionSmiles(smirks);
		} catch (InvalidSmilesException e) {
			return false;
		}
		return true;
	}
	private static boolean isKekulizationValid(IAtomContainer ac) {
		return tools.isValencyOfMoleculeValid(ac, false);
		
		
	}

}
