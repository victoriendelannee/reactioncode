/*
 * RXNFileManipulator.java
 *
 * Created on September 12, 2007, 3:42 AM
 *
 *
 * @author Syed Asad Rahman, EMBL-EBI, Cambridge, UK
 * @contact asad@ebi.ac.uk
 *
 */
package com.nih.parser;

import java.io.File;
import java.io.IOException;
import static java.lang.Integer.parseInt;
import static java.lang.System.err;
import static java.lang.System.out;
import java.util.BitSet;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;
import static java.util.logging.Level.SEVERE;
import java.util.logging.Logger;
import static java.util.logging.Logger.getLogger;
import java.util.regex.Pattern;
import static java.util.regex.Pattern.compile;
import org.openscience.cdk.AtomContainer;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IReaction;
import static org.openscience.cdk.interfaces.IReaction.Direction.BACKWARD;
import static org.openscience.cdk.interfaces.IReaction.Direction.BIDIRECTIONAL;
import static org.openscience.cdk.interfaces.IReaction.Direction.FORWARD;
import uk.ac.ebi.reactionblast.tools.BasicDebugger;
import static uk.ac.ebi.reactionblast.tools.ExtAtomContainerManipulator.aromatizeMolecule;
import static uk.ac.ebi.reactionblast.tools.ExtAtomContainerManipulator.percieveAtomTypesAndConfigureAtoms;
import static uk.ac.ebi.reactionblast.tools.ExtAtomContainerManipulator.removeHydrogensExceptSingleAndPreserveAtomID;

/**
 *
 * @author Syed Asad Rahman, EMBL-EBI, Cambridge, UK
 * @contact asad@ebi.ac.uk
 */
public class RXNFileManipulator extends BasicDebugger {

    private static final Logger LOG = getLogger(RXNFileManipulator.class.getName());
    private Integer moleculeCounter = 0;

    /**
     *
     *
     */
    public RXNFileManipulator() {
        moleculeCounter = 0;
    }

    /**
     *
     * @param rxnFile
     * @return
     * @throws Exception
     */
    public IReaction process(File rxnFile) throws Exception {

        IReaction cleanedReaction = DefaultChemObjectBuilder.getInstance().newInstance(IReaction.class);
        Map<String, Double> _Stoichiometry = new HashMap<>();
        Set<IAtomContainer> _Metabolites = new HashSet<>();
        String FileName = rxnFile.getName();
        try {
            if (FileName.endsWith(".rxn")) {
                out.println("Reading the input RXN File");
                String RXNFileName = rxnFile.getCanonicalPath();

                RXNFileImporter rxn = new RXNFileImporter();
                rxn.readFile(RXNFileName);

                String REGEX = "\\.";
                String ReactionID = "";
                Pattern p = compile(REGEX);
                String[] fields = p.split(FileName);

                ReactionID += fields[0];

//                System.out.println("****************************");
//                System.out.println("Processing: " + ReactionID);
//                System.out.println("****************************");
                IAtomContainerSet _imoledu = rxn.getReactants();
                IAtomContainerSet _imolpro = rxn.getProducts();

                out.println("Number of Reactant " + _imoledu.getAtomContainerCount());
                out.println("Number of Product " + _imolpro.getAtomContainerCount());
                cleanedReaction.setID(ReactionID);
                _Stoichiometry.clear();

                Double _TempStoic = 0.0;

                for (int i = 0; i < _imoledu.getAtomContainerCount(); i++) {
                    IAtomContainer mol = _imoledu.getAtomContainer(i);
                    _TempStoic = 0.0;
                    IAtomContainer ac = setProperty(mol);
                    String ID = ac.getID();
                    mol = mol.getBuilder().newInstance(IAtomContainer.class, ac);
                    mol.setID(ID);
                    if (_Stoichiometry.containsKey(mol.getID())) {

                        _TempStoic = _Stoichiometry.get(mol.getID()) + 1;
                        _Stoichiometry.put(mol.getID(), _TempStoic);
                    } else {
                        _Stoichiometry.put(mol.getID(), 1.0);
                        _Metabolites.add(mol);
                    }

                }

                setReactantMolecule(cleanedReaction, _Stoichiometry, _Metabolites);

                //System.out.println();
                //System.out.println("****************************");
                //System.out.println();
                for (int i = 0; i < _imolpro.getAtomContainerCount(); i++) {
                    IAtomContainer mol = _imolpro.getAtomContainer(i);
                    IAtomContainer ac = setProperty(mol);
                    String ID = ac.getID();
                    mol = mol.getBuilder().newInstance(IAtomContainer.class, ac);
                    mol.setID(ID);
                    _TempStoic = 0.0;

                    if (_Stoichiometry.containsKey(mol.getID())) {
                        _TempStoic = _Stoichiometry.get(mol.getID()) + 1;
                        _Stoichiometry.put(mol.getID(), _TempStoic);
                    } else {
                        _Stoichiometry.put(mol.getID(), 1.0);
                        _Metabolites.add(mol);
                    }

                }

                setProductMolecule(cleanedReaction, _Stoichiometry, _Metabolites);
                //As per CDK BIDIRECTION 1, Forward 2, Backward 0
                cleanedReaction.setDirection(BIDIRECTIONAL);
            }//end of if

        } catch (IOException ex) {
            err.println("Error: Generation of FingerPrint Failed the AtomContainer: ");
            ex.printStackTrace();
        }

        out.println("Reaction is read and initialized");
        return cleanedReaction;

    }

    /**
     *
     * @param rxnFile
     * @return
     * @throws Exception
     */
    public IReaction processMACiE(File rxnFile) throws Exception {

        IReaction IR = DefaultChemObjectBuilder.getInstance().newInstance(IReaction.class);
        Map<String, Double> _Stoichiometry = new HashMap<>();
        Set<IAtomContainer> _Metabolites = new HashSet<>();
        String FileName = rxnFile.getName();
        if (FileName.endsWith(".rxn")) {
            try {

                String RXNFileName = rxnFile.getCanonicalPath();
//                    molID = null;
                //System.out.println();
                RXNFileImporter rxn = new RXNFileImporter();
                rxn.readFile(RXNFileName);
                String REGEX = "\\.";
                String ReactionID = "R";
                Pattern p = compile(REGEX);
                String[] fields = p.split(FileName);
//                System.out.println(fields);
                if (!FileName.endsWith(".ov.rxn")) {
                    out.println("Reading the input MACiE rxn (stage reaction) File");
                    ReactionID = ReactionID + fields[0] + fields[1].replace("stg", ".");
                } else {
                    out.println("Reading the input MACiE rxn (overall reaction) File");
                    ReactionID = fields[0] + fields[1].replace("ov", "");
                }
//                System.out.println("****************************");
//                System.out.println("Processing: " + ReactionID);
//                System.out.println("****************************");

                IAtomContainerSet _imoledu = rxn.getReactants();
                IAtomContainerSet _imolpro = rxn.getProducts();
//
//                System.out.println(_imoledu.getAtomContainerCount());
//                System.out.println(_imolpro.getAtomContainerCount());

                IR.setID(ReactionID);

                Double _TempStoic = 0.0;
                for (int i = 0; i < _imoledu.getAtomContainerCount(); i++) {
                    IAtomContainer mol = _imoledu.getAtomContainer(i);
                    _TempStoic = 0.0;
                    IAtomContainer ac = setProperty(mol);
                    String ID = ac.getID();
                    mol = mol.getBuilder().newInstance(IAtomContainer.class, ac);
                    mol.setID(ID);
                    if (_Stoichiometry.containsKey(mol.getID())) {
                        _TempStoic = _Stoichiometry.get(mol.getID()) + 1;
                        _Stoichiometry.put(mol.getID(), _TempStoic);
                    } else {
                        _Stoichiometry.put(mol.getID(), 1.0);
                        _Metabolites.add(mol);
                    }
                }
//                System.out.println(_Metabolites.size());
                setReactantMolecule(IR, _Stoichiometry, _Metabolites);
//                System.out.println();
//                System.out.println("****************************");
//                System.out.println();
                for (int i = 0; i < _imolpro.getAtomContainerCount(); i++) {
                    IAtomContainer mol = _imolpro.getAtomContainer(i);
                    IAtomContainer ac = setProperty(mol);
                    String ID = ac.getID();
                    mol = mol.getBuilder().newInstance(IAtomContainer.class, ac);
                    mol.setID(ID);
                    _TempStoic = 0.0;
                    if (_Stoichiometry.containsKey(mol.getID())) {
                        _TempStoic = _Stoichiometry.get(mol.getID()) + 1;
                        _Stoichiometry.put(mol.getID(), _TempStoic);
                    } else {
                        _Stoichiometry.put(mol.getID(), 1.0);
                        _Metabolites.add(mol);
                    }
                }
                setProductMolecule(IR, _Stoichiometry, _Metabolites);
                //As per CDK BIDIRECTION 1, Forward 2, Backward 0
                IR.setDirection(BIDIRECTIONAL);
                //Print the reaction
//                printReaction(cleanedReaction);
            } catch (IOException ex) {
                getLogger(RXNFileManipulator.class.getName()).log(SEVERE, null, ex);
            }
        }//end of if

        out.println("Reaction is read and initialized");
        return IR;
    }

    /**
     *
     * @param rxnFile
     * @return
     * @throws Exception
     */
    public IReaction processIntEnz(File rxnFile) throws Exception {

        IReaction IR = DefaultChemObjectBuilder.getInstance().newInstance(IReaction.class);
        Map<String, Double> _Stoichiometry = new HashMap<>();
        Set<IAtomContainer> _Metabolites = new HashSet<>();
        String FileName = rxnFile.getName();
        try {
            if (FileName.endsWith(".rxn")) {
                out.println("Reading the input RehA rxn File");
                String RXNFileName = rxnFile.getCanonicalPath();
                RXNFileImporter rxn = new RXNFileImporter();
                rxn.readFile(RXNFileName);
                String REGEX = "\\.";
                String ReactionID = "";
                Pattern p = compile(REGEX);
                String[] fields = p.split(FileName);
                ReactionID += fields[0];
                int direction = parseInt(fields[1]);
//                System.out.println("****************************");
//                System.out.println("Processing: " + ReactionID);
//                System.out.println("****************************");

                IAtomContainerSet _imoledu = rxn.getReactants();
                IAtomContainerSet _imolpro = rxn.getProducts();
                IR.setID(ReactionID);
                _Stoichiometry.clear();
                Double _TempStoic = 0.0;
                for (int i = 0; i < _imoledu.getAtomContainerCount(); i++) {
                    IAtomContainer mol = _imoledu.getAtomContainer(i);
                    _TempStoic = 0.0;
                    IAtomContainer ac = setProperty(mol);
                    String ID = ac.getID();
                    mol = mol.getBuilder().newInstance(IAtomContainer.class, ac);
                    mol.setID(ID);
                    if (_Stoichiometry.containsKey(mol.getID())) {
                        _TempStoic = _Stoichiometry.get(mol.getID()) + 1;
                        _Stoichiometry.put(mol.getID(), _TempStoic);
                    } else {
                        _Stoichiometry.put(mol.getID(), 1.0);
                        _Metabolites.add(mol);
                    }
                }
                setReactantMolecule(IR, _Stoichiometry, _Metabolites);
                //System.out.println();
                //System.out.println("****************************");
                //System.out.println();
                for (int i = 0; i < _imolpro.getAtomContainerCount(); i++) {
                    IAtomContainer mol = _imolpro.getAtomContainer(i);
                    IAtomContainer ac = setProperty(mol);
                    String ID = ac.getID();
                    mol = mol.getBuilder().newInstance(IAtomContainer.class, ac);
                    mol.setID(ID);
                    _TempStoic = 0.0;
                    if (_Stoichiometry.containsKey(mol.getID())) {
                        _TempStoic = _Stoichiometry.get(mol.getID()) + 1;
                        _Stoichiometry.put(mol.getID(), _TempStoic);
                    } else {
                        _Stoichiometry.put(mol.getID(), 1.0);
                        _Metabolites.add(mol);
                    }
                }
                setProductMolecule(IR, _Stoichiometry, _Metabolites);
                //As per IntEnz 0 for undefined direction, 1 for LR, 2 for RL and 3 for bidirectional
                //As per CDK BIDIRECTION 1, Forward 2, Backward 0
                if (direction == 1) {
                    IR.setDirection(FORWARD);
                } else if (direction == 2) {
                    IR.setDirection(BACKWARD);
                } else if (direction == 3 || direction == 0) {
                    IR.setDirection(BIDIRECTIONAL);
                }

            } //end of if
        } catch (IOException ex) {
            getLogger(RXNFileManipulator.class.getName()).log(SEVERE, null, ex);
        }
        out.println("Reaction is read and initialized");
        return IR;
//
    }

    private void setReactantMolecule(IReaction IR, Map<String, Double> _Stoichiometry, Set<IAtomContainer> _Metabolites) {
////
//        System.out.println("-------------------");
//        System.out.println("Reactants");
//        System.out.println("-------------------");

        Iterator<IAtomContainer> it = _Metabolites.iterator();
//        System.out.println("Stoic Map Size: " + _Stoichiometry.size());
        while (it.hasNext()) {
            IAtomContainer mol = it.next();
//            System.out.println("Insertion: " + mol.getID() + ", " + _Stoichiometry.get(mol.getID())
//                    + ", ac:" + mol.getAtomCount());
//            mol.setProperty("STOICHIOMETRY", _Stoichiometry.get(mol.getID()));
            IR.addReactant(mol, _Stoichiometry.get(mol.getID()));
        }

        _Metabolites.clear();
        _Stoichiometry.clear();
    }

    private void setProductMolecule(IReaction IR, Map<String, Double> _Stoichiometry, Set<IAtomContainer> _Metabolites) {

//        System.out.println("-------------------");
//        System.out.println("Products");
//        System.out.println("-------------------");
        Iterator<IAtomContainer> it = _Metabolites.iterator();

        while (it.hasNext()) {

            IAtomContainer mol = it.next();
//            System.out.println("Insertion: " + mol.getID() + ", " + _Stoichiometry.get(mol.getID())
//                    + ", ac:" + mol.getAtomCount());
//            mol.setProperty("STOICHIOMETRY", _Stoichiometry.get(mol.getID()));
            IR.addProduct(mol, _Stoichiometry.get(mol.getID()));
        }

        _Metabolites.clear();
        _Stoichiometry.clear();
    }

    private IAtomContainer setProperty(IAtomContainer mol) throws Exception {
        String molID = mol.getID();
//         System.out.println("After Cleanup");
//        System.out.println(CDKChemaxonIOConveter.getChemAxonMolecule(mol).exportToFormat("mol:V2"));

        IAtomContainer molecule = new AtomContainer(mol);
        removeHydrogensExceptSingleAndPreserveAtomID(mol);
        percieveAtomTypesAndConfigureAtoms(molecule);
        aromatizeMolecule(molecule);
        molecule.setID(molID);

       
//        printAtoms(molecule);
//        System.out.println(molecule.getID());
        return molecule;
    }
}
