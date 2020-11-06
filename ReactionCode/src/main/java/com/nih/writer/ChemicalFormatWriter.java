package com.nih.writer;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IReaction;
import org.openscience.cdk.DefaultChemObjectBuilder;

import com.nih.tools.tools;

public class ChemicalFormatWriter {

	/**
	 * Write RDF
	 * @param reactions (/home/mydirectory/myreaction.rdf)
	 * @param path
	 * @throws IOException 
	 */
	public void writeRDF(IReaction reactions, String path) throws IOException {
		File file = new File(path);
		try (MDLV2000RDFWriter writer = new MDLV2000RDFWriter(new FileWriter(file))) {
			writer.write(reactions);
			//writer.setRdFieldsProperties(map);
			writer.close();
		} catch (CDKException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	/**
	 * Write RXN
	 * @param reactions (/home/mydirectory/myreaction.rxn)
	 * @param path
	 * @throws IOException 
	 */
	public void writeRXN(IReaction reactions, String path) throws IOException {
		File file = new File(path);
		try (MDLV2000RXNWriter writer = new MDLV2000RXNWriter(new FileWriter(file))) {
			writer.write(reactions);
			//writer.setRdFieldsProperties(map);
			writer.close();
		} catch (CDKException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	/**
	 * Write SMIRKS (=Reaction SMARTS)
	 * @param reaction
	 * @return
	 * @throws CDKException 
	 */
	public String writeSMIRKS(IReaction reaction) throws CDKException {
		return tools.makeSmiles(reaction, true);
	}
	
	/**
	 * Write Reaction SMILES
	 * @param reaction
	 * @return
	 * @throws CDKException 
	 */
	public String writeReactionSMILES(IReaction reaction) throws CDKException {
		return tools.makeSmiles(reaction, false);
	}
	
	/**
	 * Write SMARTS
	 * @param reaction
	 * @return
	 * @throws CDKException 
	 */
	public String writeSMARTS(IAtomContainer ac) throws CDKException {
		return tools.makeSmiles(ac, true, true, false);
	}
	
	/**
	 * Write SMARTS
	 * @param reaction
	 * @return
	 * @throws CDKException 
	 */
	public String writeSMARTS(IAtomContainerSet set) throws CDKException {
    	IAtomContainer aggregate = DefaultChemObjectBuilder.getInstance().newAtomContainer();

    	for (IAtomContainer ac : set.atomContainers()) {
    		aggregate.add(ac);
    	}
		return tools.makeSmiles(aggregate, true, true, false);
	}
	
	/**
	 * Write SMILES
	 * @param reaction
	 * @return
	 * @throws CDKException 
	 */
	public String writeSMILES(IAtomContainer ac) throws CDKException {
		return tools.makeSmiles(ac, true, false, false);
	}
	
	/**
	 * Write SMILES
	 * @param reaction
	 * @return
	 * @throws CDKException 
	 */
	public String writeSMILES(IAtomContainerSet set) throws CDKException {
    	IAtomContainer aggregate = DefaultChemObjectBuilder.getInstance().newAtomContainer();

    	for (IAtomContainer ac : set.atomContainers()) {
    		aggregate.add(ac);
    	}
		return tools.makeSmiles(aggregate, true, false, false);
	}
}
