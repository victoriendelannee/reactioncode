package com.nih.parser;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.Reader;
import java.io.StringReader;

import static java.lang.Integer.valueOf;
import static java.lang.System.getProperty;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.StringTokenizer;
//import java.util.logging.Logger;
//import static java.util.logging.Logger.getLogger;
import static org.openscience.cdk.CDKConstants.ATOM_ATOM_MAPPING;
import static org.openscience.cdk.CDKConstants.TITLE;

import org.openscience.cdk.aromaticity.Aromaticity;
import org.openscience.cdk.aromaticity.ElectronDonation;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.graph.CycleFinder;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IChemFile;
import org.openscience.cdk.interfaces.IChemModel;
import org.openscience.cdk.interfaces.IChemObject;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.interfaces.IChemSequence;
import org.openscience.cdk.interfaces.IMapping;
import org.openscience.cdk.interfaces.IReaction;
import org.openscience.cdk.interfaces.IReactionSet;
import org.openscience.cdk.io.DefaultChemObjectReader;
import static org.openscience.cdk.io.IChemObjectReader.Mode.RELAXED;
import org.openscience.cdk.io.formats.IResourceFormat;
import static org.openscience.cdk.io.formats.MDLRXNFormat.getInstance;
import org.openscience.cdk.tools.ILoggingTool;

import com.nih.tools.ColouredSystemOutPrintln;
import com.nih.writer.MDLV2000RDFWriter;

import static org.openscience.cdk.tools.LoggingToolFactory.createLoggingTool;

/**
 * Reads a molecule from an MDL RXN and MDLRD file {
 *
 * @cdk.cite DAL92}. This MDL RXN reader uses the MDLV2000 reader to read each
 * mol file
 * @cdk.module io
 * 
 *
 * @author Egon Willighagen
 * @author Thomas Kuhn
 * @cdk.created 2003-07-24
 *
 * @cdk.keyword file format: MDL RXN and MDL RDF
 * @cdk.bug 1849923
 */
public class MDLReactionsFileReader extends DefaultChemObjectReader {

	private boolean aromatize = false;
	private boolean split = false;
	private String filename;
	int reactionCount = 0;
	List<File> files = new ArrayList<File>();
	
	private static ILoggingTool logger
	= createLoggingTool(MDLReactionsFileReader.class);
	//private static Logger LOG = getLogger(MDLReactionsFileReader.class.getName());
	BufferedReader input = null;

	/**
	 * Constructs a new MDLReader that can read AtomContainer from a given
	 * Reader.
	 *
	 * @param in The Reader to read from
	 */
	public MDLReactionsFileReader(Reader in) {
		this(in, RELAXED);
	}

	/**
	 *
	 * @param in
	 * @param mode
	 */
	public MDLReactionsFileReader(Reader in, Mode mode) {
		if (in instanceof BufferedReader) {
			input = (BufferedReader) in;
		} else {
			input = new BufferedReader(in);
		}
		super.mode = mode;
	}

	/**
	 *
	 * @param input
	 */
	public MDLReactionsFileReader(InputStream input) {
		this(input, RELAXED);
	}

	/**
	 *
	 * @param input
	 * @param mode
	 */
	public MDLReactionsFileReader(InputStream input, Mode mode) {
		this(new InputStreamReader(input), mode);
	}

	/**
	 *
	 */
	public MDLReactionsFileReader() {
		this(new StringReader(""));
	}

	/**
	 *
	 * @return
	 */
	@Override
	public IResourceFormat getFormat() {
		return getInstance();
	}

	/**
	 *
	 * @param input
	 * @throws CDKException
	 */
	@Override
	public void setReader(Reader input) throws CDKException {
		if (input instanceof BufferedReader) {
			this.input = (BufferedReader) input;
		} else {
			this.input = new BufferedReader(input);
		}
	}

	/**
	 *
	 * @param input
	 * @throws CDKException
	 */
	@Override
	public void setReader(InputStream input) throws CDKException {
		setReader(new InputStreamReader(input));
	}

	/**
	 *
	 * @param classObject
	 * @return
	 */
	@Override
	public boolean accepts(Class classObject) {
		Class[] interfaces = classObject.getInterfaces();
		for (Class intf : interfaces) {
			if (IChemModel.class.equals(intf)) {
				return true;
			}
			if (IChemFile.class.equals(intf)) {
				return true;
			}
			if (IReaction.class.equals(intf)) {
				return true;
			}
		}
		Class superClass = classObject.getSuperclass();
		if (superClass != null) {
			return this.accepts(superClass);
		}
		return false;
	}

	//    /**
	//     * Takes an object which subclasses IChemObject, e.g.AtomContainer, and will
	//     * read this (from file, database, internet etc). If the specific
	//     * implementation does not support a specific IChemObject it will throw an
	//     * Exception.
	//     *
	//     * @param <T>
	//     * @param object The object that subclasses IChemObject
	//     * @return The IChemObject read
	//     * @exception CDKException
	//     */
	//    @Override
	//    public <T extends IChemObject> T read(T object) throws CDKException {
	//        if (object instanceof IReaction) {
	//            return (T) readReaction(object.getBuilder());
	//        } else if (object instanceof IReactionSet) {
	//            IReactionSet reactionSet = object.getBuilder().newInstance(IReactionSet.class);
	//            reactionSet.addReaction(readReaction(object.getBuilder()));
	//            return (T) reactionSet;
	//        } else if (object instanceof IChemModel) {
	//            IChemModel model = object.getBuilder().newInstance(IChemModel.class);
	//            IReactionSet reactionSet = object.getBuilder().newInstance(IReactionSet.class);
	//            reactionSet.addReaction(readReaction(object.getBuilder()));
	//            model.setReactionSet(reactionSet);
	//            return (T) model;
	//        } else if (object instanceof IChemFile) {
	//            IChemFile chemFile = object.getBuilder().newInstance(IChemFile.class);
	//            IChemSequence sequence = object.getBuilder().newInstance(IChemSequence.class);
	//            sequence.addChemModel(read(object.getBuilder().newInstance(IChemModel.class)));
	//            chemFile.addChemSequence(sequence);
	//            return (T) chemFile;
	//        } else {
	//            throw new CDKException("Only supported are Reaction and ChemModel, and not "
	//                    + object.getClass().getName() + ".");
	//        }
	//    }
	/**
	 * Takes an object which subclasses IChemObject, e.g.AtomContainer, and will
	 * read this (from file, database, internet etc). If the specific
	 * implementation does not support a specific IChemObject it will throw an
	 * Exception.
	 *
	 * @param <T>
	 * @param object The object that subclasses IChemObject
	 * @return The IChemObject read
	 * @exception CDKException
	 */
	@Override
	public <T extends IChemObject> T read(T object) throws CDKException {
		if (object instanceof IReaction) {
			return (T) readReaction(object.getBuilder());
		} else if (object instanceof IReactionSet) {
			IReactionSet reactionSet = object.getBuilder().newInstance(IReactionSet.class);
			return (T) readReactions(object.getBuilder());
		} else if (object instanceof IChemModel) {
			IChemModel model = object.getBuilder().newInstance(IChemModel.class);
			IReactionSet reactionSet = object.getBuilder().newInstance(IReactionSet.class);
			model.setReactionSet((IReactionSet) readReactions(object.getBuilder()));
			return (T) model;
		} else if (object instanceof IChemFile) {
			IChemFile chemFile = object.getBuilder().newInstance(IChemFile.class);
			IChemSequence sequence = object.getBuilder().newInstance(IChemSequence.class);
			sequence.addChemModel(read(object.getBuilder().newInstance(IChemModel.class)));
			chemFile.addChemSequence(sequence);
			return (T) chemFile;
		} else {
			throw new CDKException("Only supported are Reaction and ChemModel, and not "
					+ object.getClass().getName() + ".");
		}
	}
	/**
	 *
	 * @param object
	 * @return
	 */
	public boolean accepts(IChemObject object) {
		if (object instanceof IReaction) {
			return true;
		} else if (object instanceof IChemModel) {
			return true;
		} else if (object instanceof IChemFile) {
			return true;
		} else if (object instanceof IReactionSet) {
			return true;
		}
		return false;
	}

	/**
	 * Read a Reaction from a file in MDL RXN format
	 *
	 * @return The Reaction that was read from the MDL file.
	 */
	private IReaction readReaction2(IChemObjectBuilder builder) throws CDKException {
		IReaction reaction = builder.newInstance(IReaction.class);
		//if aromatize
		ElectronDonation model       = ElectronDonation.daylight();
		CycleFinder      cycles      = Cycles.or(Cycles.all(), Cycles.all(6));
		Aromaticity      aromaticity = new Aromaticity(model, cycles);
		String countsLine = null;
		int count = 0;
		try {
			String line = input.readLine();
			if (line == null)
				return null;
			while (!line.contains("$RXN")) {
				line = input.readLine();
				if (line.contains("RIREG")) {
					String[] t = line.split("$RIREG");
					for (int i = t.length - 1; i > -1; i--) {
						String s = t[i].replace(" ", "");
						if (! s.equals("")) {
							reaction.setID(s);
							break;
						}
					}
				}
			}
			while (count != 2) {
				line = input.readLine();
				if (line.contains("[a-zA-Z]+") == false) {
					List<String> t = Arrays.asList(line.split(" "));
					for (int i = 0; i < t.size(); i ++) {
						if (t.get(i).matches("-?\\d+(\\.\\d+)?")) {
							count ++;
						}
					}
					if (count == 2) countsLine = line;
				}
			}

			//            input.readLine(); // first line should be $RXN
			//            input.readLine(); // second line
			//            input.readLine(); // third line
			//            input.readLine(); // fourth line
		} catch (IOException exception) {
			logger.debug(exception);
			throw new CDKException("Error while reading header of RXN file", exception);
		}

		int reactantCount = 0;
		int productCount = 0;
		try {
			//String countsLine = input.readLine();
			/* this line contains the number of reactants
             and products */
			StringTokenizer tokenizer = new StringTokenizer(countsLine);
			reactantCount = valueOf(tokenizer.nextToken());
			logger.info("Expecting " + reactantCount + " reactants in file");
			productCount = valueOf(tokenizer.nextToken());
			logger.info("Expecting " + productCount + " products in file");
		} catch (NumberFormatException exception) {
			logger.debug(exception);
			throw new CDKException("Error while counts line of RXN file", exception);
		}

		// now read the reactants
		try {
			for (int i = 1; i <= reactantCount; i++) {
				StringBuilder molFile = new StringBuilder();
				String line = input.readLine(); // announceMDLFileLine
				String molFileLine = "";
				if (line.contains("$MOL")) {
					do {
						molFileLine = input.readLine();
						molFile.append(molFileLine);
						molFile.append(getProperty("line.separator"));
					} while (!molFileLine.equals("M  END"));
				}

				// read MDL molfile content
				// Changed this to mdlv2000 reader
				//                MDLV2000Reader reader = new MDLV2000Reader(
				//                        new StringReader(molFile.toString()),
				//                        super.mode);
				//                IAtomContainer reactant = reader.read(
				//                        builder.newInstance(IAtomContainer.class));
				MDLFileReader reader = new MDLFileReader(new StringReader(molFile.toString()),
						super.mode);
				IAtomContainer reactant = reader.getAtomContainer();
				if (reactant == null) {
					continue;
				}
				if (aromatize) {
					 aromaticity.apply(reactant);
				}
				// add reactant mol ID
				String readMolID = (String) reactant.getProperty(TITLE);
				if (readMolID != null) {
					reactant.setID(readMolID.trim());
				}
				// add reactant
				reaction.addReactant(reactant);
			}
		} 
		//        catch (CDKException exception) {
		//            // rethrow exception from MDLReader
		//            throw exception;
		//        } 
		catch (IOException | IllegalArgumentException exception) {
			logger.debug(exception);
			throw new CDKException("Error while reading reactant", exception);
		}

		// now read the products
		try {
			for (int i = 1; i <= productCount; i++) {
				StringBuilder molFile = new StringBuilder();
				String line = input.readLine(); // String announceMDLFileLine = 
				String molFileLine = "";
				if (line.contains("$MOL")) {
					do {
						molFileLine = input.readLine();
						molFile.append(molFileLine);
						molFile.append(getProperty("line.separator"));
					} while (!molFileLine.equals("M  END"));
				}

				// read MDL molfile content
				//                MDLV2000Reader reader = new MDLV2000Reader(
				//                        new StringReader(molFile.toString()));
				//                IAtomContainer product = reader.read(
				//                        builder.newInstance(IAtomContainer.class));
				MDLFileReader reader = new MDLFileReader(new StringReader(molFile.toString()),
						super.mode);
				IAtomContainer product = reader.getAtomContainer();

				if (product == null) {
					continue;
				}
				if (aromatize) {
					 aromaticity.apply(product);
				}

				// add product molID
				String readMolID = (String) product.getProperty(TITLE);
				if (readMolID != null) {
					product.setID(readMolID.trim());
				}
				// add product
				reaction.addProduct(product);
			}
		} 
		//        catch (CDKException exception) {
		//            // rethrow exception from MDLReader
		//            throw exception;
		//        } 
		catch (IOException | IllegalArgumentException exception) {
			logger.debug(exception);
			throw new CDKException("Error while reading products", exception);
		}

		// now try to map things, if wanted
		logger.info("Reading atom-atom mapping from file");
		// distribute all atoms over two GraphAtomContainer's
		IAtomContainer reactingSide = builder.newInstance(IAtomContainer.class);
		java.util.Iterator<IAtomContainer> molecules = reaction.getReactants().atomContainers().iterator();
		while (molecules.hasNext()) {
			reactingSide.add(molecules.next());
		}
		IAtomContainer producedSide = builder.newInstance(IAtomContainer.class);
		molecules = reaction.getProducts().atomContainers().iterator();
		while (molecules.hasNext()) {
			producedSide.add(molecules.next());
		}

		// map the atoms
		int mappingCount = 0;
		for (int i = 0; i < reactingSide.getAtomCount(); i++) {
			for (int j = 0; j < producedSide.getAtomCount(); j++) {
				IAtom eductAtom = reactingSide.getAtom(i);
				IAtom productAtom = producedSide.getAtom(j);
				if (eductAtom.getProperty(ATOM_ATOM_MAPPING) != null
						&& eductAtom.getProperty(ATOM_ATOM_MAPPING).equals(productAtom.getProperty(ATOM_ATOM_MAPPING))) {
					reaction.addMapping(
							builder.newInstance(IMapping.class, eductAtom, productAtom));
					mappingCount++;
					break;
				}
			}
		}
		logger.info("Mapped atom pairs: " + mappingCount);

		return reaction;
	}

	private IReaction readReaction(IChemObjectBuilder builder) throws CDKException {
		IReaction reaction = builder.newInstance(IReaction.class);
		String line = "";
		int reactantCount = 0;
		int productCount = 0;
		
		//if aromatize
		ElectronDonation model       = ElectronDonation.daylight();
		CycleFinder      cycles      = Cycles.or(Cycles.all(), Cycles.all(6));
		Aromaticity      aromaticity = new Aromaticity(model, cycles);
		
		try {
			while ((line = input.readLine()) != null){

				if (line == null &&  reaction.getReactantCount() == reactantCount && reaction.getProductCount() != productCount)
					return reaction;
				else if (line == null && (reaction.getReactantCount() != reactantCount || reaction.getProductCount() != productCount))
					return null;
				
				while (!line.contains("$RXN")) 
					line = input.readLine();
				
				
					int count = 0;
					String countsLine = null;
					int headerCount = 1;
					try {
//						while (!line.contains("$RXN")) {
//							line = input.readLine();
//							if (line == null) return reactionSet;
//						}
						while (count != 2) {
							line = input.readLine();
							headerCount++;
							if (line.contains("[a-zA-Z]+") == false && headerCount > 3) {
								List<String> t = Arrays.asList(line.split(" "));
								for (int i = 0; i < t.size(); i ++) {
									if (t.get(i).matches("-?\\d+(\\.\\d+)?")) {
										count++;
									}
								}
								if (count == 2) countsLine = line;
								else count = 0;
							}
							if (count != 2 && line.matches("^[a-zA-Z0-9]*$") && headerCount < 4)
								reaction.setID(line);
						}

						//                input.readLine(); // first line should be $RXN
						//                input.readLine(); // second line
						//                input.readLine(); // third line
						//                input.readLine(); // fourth line
					} catch (IOException exception) {
						logger.debug(exception);
						throw new CDKException("Error while reading header of RXN file", exception);
					}

					try {

						//String countsLine = input.readLine();
						/* this line contains the number of reactants
					         and products */
						StringTokenizer tokenizer = new StringTokenizer(countsLine);
						reactantCount = valueOf(tokenizer.nextToken());
						logger.info("Expecting " + reactantCount + " reactants in file");
						productCount = valueOf(tokenizer.nextToken());
						logger.info("Expecting " + productCount + " products in file");
					} catch (NumberFormatException exception) {
						logger.debug(exception);
						throw new CDKException("Error while counts line of RXN file", exception);
					}

					// now read the reactants
					try {
						for (int i = 1; i <= reactantCount; i++) {
							StringBuilder molFile = new StringBuilder();
							line = input.readLine(); // announceMDLFileLine
							String molFileLine = "";
							if (line.contains("$MOL")) {
								do {
									molFileLine = input.readLine();
									molFile.append(molFileLine);
									molFile.append(getProperty("line.separator"));
								} while (!molFileLine.equals("M  END"));
							}

							// read MDL molfile content
							// Changed this to mdlv2000 reader
							//                    MDLV2000Reader reader = new MDLV2000Reader(
							//                            new StringReader(molFile.toString()),
							//                            super.mode);
							//                    IAtomContainer reactant = reader.read(
							//                            builder.newInstance(IAtomContainer.class));
							MDLFileReader reader = new MDLFileReader(new StringReader(molFile.toString()),
									super.mode);
							IAtomContainer reactant = reader.getAtomContainer();
							if (reactant.getProperty("agent").equals(true)) {
								reaction.setProperty("agent", true);
							}
							if (reactant == null) {
								continue;
							}
							if (aromatize) {
								 aromaticity.apply(reactant);
							}
							// add reactant mol ID
							String readMolID = (String) reactant.getProperty(TITLE);
							if (readMolID != null) {
								reactant.setID(readMolID.trim());
							}
							// add reactant
							reaction.addReactant(reactant);
						}
					} 
					//            catch (CDKException exception) {
					//                // rethrow exception from MDLReader
					//                throw exception;
					//            } 
					catch (IOException | IllegalArgumentException exception) {
						logger.debug(exception);
						throw new CDKException("Error while reading reactant", exception);
					}

					// now read the products
					try {
						for (int i = 1; i <= productCount; i++) {
							StringBuilder molFile = new StringBuilder();
							line = input.readLine(); // String announceMDLFileLine = 
							String molFileLine = "";
							if (line.contains("$MOL")) {
								do {
									molFileLine = input.readLine();
									molFile.append(molFileLine);
									molFile.append(getProperty("line.separator"));
								} while (!molFileLine.equals("M  END"));
							}

							// read MDL molfile content
							//                    MDLV2000Reader reader = new MDLV2000Reader(
							//                            new StringReader(molFile.toString()));
							//                    IAtomContainer product = reader.read(
							//                            builder.newInstance(IAtomContainer.class));
							MDLFileReader reader = new MDLFileReader(new StringReader(molFile.toString()),
									super.mode);
							IAtomContainer product = reader.getAtomContainer();

							if (product == null) {
								continue;
							}
							if (aromatize) {
								 aromaticity.apply(product);
							}

							// add product molID
							String readMolID = (String) product.getProperty(TITLE);
							if (readMolID != null) {
								product.setID(readMolID.trim());
							}
							// add product
							reaction.addProduct(product);
						}
					} 
					//            catch (CDKException exception) {
					//                // rethrow exception from MDLReader
					//                throw exception;
					//            } 
					catch (IOException | IllegalArgumentException exception) {
						logger.debug(exception);
						throw new CDKException("Error while reading products", exception);
					}

					if (reaction.getProperty("agent") == null) {
						reaction.setProperty("agent", false);
					}

					// now try to map things, if wanted
					logger.info("Reading atom-atom mapping from file");
					// distribute all atoms over two GraphAtomContainer's
					IAtomContainer reactingSide = builder.newInstance(IAtomContainer.class);
					java.util.Iterator<IAtomContainer> molecules = reaction.getReactants().atomContainers().iterator();
					while (molecules.hasNext()) {
						reactingSide.add(molecules.next());
					}
					IAtomContainer producedSide = builder.newInstance(IAtomContainer.class);
					molecules = reaction.getProducts().atomContainers().iterator();
					while (molecules.hasNext()) {
						producedSide.add(molecules.next());
					}

					// map the atoms
					int mappingCount = 0;
					for (int i = 0; i < reactingSide.getAtomCount(); i++) {
						for (int j = 0; j < producedSide.getAtomCount(); j++) {
							IAtom eductAtom = reactingSide.getAtom(i);
							IAtom productAtom = producedSide.getAtom(j);
							if (eductAtom.getProperty(ATOM_ATOM_MAPPING) != null
									&& eductAtom.getProperty(ATOM_ATOM_MAPPING).equals(productAtom.getProperty(ATOM_ATOM_MAPPING))) {
								reaction.addMapping(
										builder.newInstance(IMapping.class, eductAtom, productAtom));
								mappingCount++;
								break;
							}
						}
					}
					logger.info("Mapped atom pairs: " + mappingCount);
					//						reactionSet.addReaction(reaction);
//				}
			}

		} catch (IllegalArgumentException | IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		return reaction;
	}
	
	/**
	 * Read a Reaction from a file in MDL RXN format
	 *
	 * @return The Reaction that was read from the MDL file.
	 */
	private IReactionSet readReactions(IChemObjectBuilder builder) throws CDKException {
		files = new ArrayList<File>();
		IReactionSet reactionSet = builder.newInstance(IReactionSet.class);
		IReaction reaction = null;
		Map<Object,Object> properties = new HashMap<Object,Object>();
		
		String line = "";
		reactionCount = 0;
		
		//if aromatize
		ElectronDonation model       = ElectronDonation.daylight();
		CycleFinder      cycles      = Cycles.or(Cycles.all(), Cycles.all(6));
		Aromaticity      aromaticity = new Aromaticity(model, cycles);
		
		try {
			while ((line = input.readLine()) != null){
				do {
					if (line.contains("$DTYPE")) {
						properties.put(line.replace("$DTYPE ", "").toString(), input.readLine().replace("$DATUM ", "").toString());
					}
					if (line.contains("$RFMT")) {
						if (reaction != null) {
							if (split == true && reactionSet.getReactionCount() == 10000) {
								writeRDFChunk(reactionSet, filename);
								reactionSet.removeAllReactions();
								System.out.println(ColouredSystemOutPrintln.ANSI_GREEN + reactionCount + 
										" reactions are parsed and a tempory RD file has been created" + 
										ColouredSystemOutPrintln.ANSI_RESET);
							}
							else if (split == false && reactionSet.getReactionCount() == 10000)  {
								System.out.println(ColouredSystemOutPrintln.ANSI_GREEN + reactionCount + 
										" reactions have been read" + 
										ColouredSystemOutPrintln.ANSI_RESET);
							}
							reaction.setProperties(properties);
							reactionSet.addReaction(reaction);
							reactionCount++;
						}
						reaction = builder.newInstance(IReaction.class);
						properties = new HashMap<Object,Object>();

						if (line.contains("$RIREG")) {
							String[] t = line.split("RIREG");
							for (int i = t.length - 1; i > -1; i--) {
								String s = t[i].replace(" ", "");
								if (!s.equals("")) {
									reaction.setID(s);
									break;
								}
								else {
									reaction.setID(reactionCount+"");
								}
							}
						}
					}
					line = input.readLine();
					if (line == null) {
						reaction.setProperties(properties);
						reactionSet.addReaction(reaction);
						return reactionSet;
					}
				}
				while (!line.contains("$RXN"));

				//line == $RXN" extract reactant and product
					int count = 0;
					String countsLine = null;
					try {
						while (count != 2) {
							line = input.readLine();
							if (line.contains("[a-zA-Z]+") == false) {
								List<String> t = Arrays.asList(line.split(" "));
								for (int i = 0; i < t.size(); i ++) {
									if (t.get(i).matches("-?\\d+(\\.\\d+)?")) {
										count ++;
									}
								}
								if (count == 2) countsLine = line;
								else count = 0;
							}
						}

						//                input.readLine(); // first line should be $RXN
						//                input.readLine(); // second line
						//                input.readLine(); // third line
						//                input.readLine(); // fourth line
					} catch (IOException exception) {
						logger.debug(exception);
						throw new CDKException("Error while reading header of RXN file", exception);
					}

					int reactantCount = 0;
					int productCount = 0;
					try {

						//String countsLine = input.readLine();
						/* this line contains the number of reactants
					         and products */
						StringTokenizer tokenizer = new StringTokenizer(countsLine);
						reactantCount = valueOf(tokenizer.nextToken());
						logger.info("Expecting " + reactantCount + " reactants in file");
						productCount = valueOf(tokenizer.nextToken());
						logger.info("Expecting " + productCount + " products in file");
					} catch (NumberFormatException exception) {
						logger.debug(exception);
						throw new CDKException("Error while counts line of RXN file", exception);
					}

					// now read the reactants
					try {
						for (int i = 1; i <= reactantCount; i++) {
							StringBuilder molFile = new StringBuilder();
							line = input.readLine(); // announceMDLFileLine
							String molFileLine = "";
							if (line.contains("$MOL")) {
								do {
									molFileLine = input.readLine();
									molFile.append(molFileLine);
									molFile.append(getProperty("line.separator"));
								} while (!molFileLine.equals("M  END"));
							}

							// read MDL molfile content
							// Changed this to mdlv2000 reader
							//                    MDLV2000Reader reader = new MDLV2000Reader(
							//                            new StringReader(molFile.toString()),
							//                            super.mode);
							//                    IAtomContainer reactant = reader.read(
							//                            builder.newInstance(IAtomContainer.class));
							MDLFileReader reader = new MDLFileReader(new StringReader(molFile.toString()),
									super.mode);
							IAtomContainer reactant = reader.getAtomContainer();
							if (reactant.getProperty("agent").equals(true)) {
								reaction.setProperty("agent", true);
							}
							if (reactant == null) {
								continue;
							}
							if (aromatize) {
								 aromaticity.apply(reactant);
							}
							// add reactant mol ID
							String readMolID = (String) reactant.getProperty(TITLE);
							if (readMolID != null) {
								reactant.setID(readMolID.trim());
							}
							// add reactant
							reaction.addReactant(reactant);
						}
					} 
					//            catch (CDKException exception) {
					//                // rethrow exception from MDLReader
					//                throw exception;
					//            } 
					catch (IOException | IllegalArgumentException exception) {
						logger.debug(exception);
						throw new CDKException("Error while reading reactant", exception);
					}

					// now read the products
					try {
						for (int i = 1; i <= productCount; i++) {
							StringBuilder molFile = new StringBuilder();
							line = input.readLine(); // String announceMDLFileLine = 
							String molFileLine = "";
							if (line.contains("$MOL")) {
								do {
									molFileLine = input.readLine();
									molFile.append(molFileLine);
									molFile.append(getProperty("line.separator"));
								} while (!molFileLine.equals("M  END"));
							}

							// read MDL molfile content
							//                    MDLV2000Reader reader = new MDLV2000Reader(
							//                            new StringReader(molFile.toString()));
							//                    IAtomContainer product = reader.read(
							//                            builder.newInstance(IAtomContainer.class));
							MDLFileReader reader = new MDLFileReader(new StringReader(molFile.toString()),
									super.mode);
							IAtomContainer product = reader.getAtomContainer();
							if (product == null) {
								continue;
							}
							if (aromatize) {
								 aromaticity.apply(product);
							}
							// add product molID
							String readMolID = (String) product.getProperty(TITLE);
							if (readMolID != null) {
								product.setID(readMolID.trim());
							}
							// add product
							reaction.addProduct(product);
						}
					} 
					//            catch (CDKException exception) {
					//                // rethrow exception from MDLReader
					//                throw exception;
					//            } 
					catch (IOException | IllegalArgumentException exception) {
						logger.debug(exception);
						throw new CDKException("Error while reading products", exception);
					}

					if (reaction.getProperty("agent") == null) {
						reaction.setProperty("agent", false);
					}

					// now try to map things, if wanted
					logger.info("Reading atom-atom mapping from file");
					// distribute all atoms over two GraphAtomContainer's
					IAtomContainer reactingSide = builder.newInstance(IAtomContainer.class);
					java.util.Iterator<IAtomContainer> molecules = reaction.getReactants().atomContainers().iterator();
					while (molecules.hasNext()) {
						reactingSide.add(molecules.next());
					}
					IAtomContainer producedSide = builder.newInstance(IAtomContainer.class);
					molecules = reaction.getProducts().atomContainers().iterator();
					while (molecules.hasNext()) {
						producedSide.add(molecules.next());
					}

					// map the atoms
					int mappingCount = 0;
					for (int i = 0; i < reactingSide.getAtomCount(); i++) {
						for (int j = 0; j < producedSide.getAtomCount(); j++) {
							IAtom eductAtom = reactingSide.getAtom(i);
							IAtom productAtom = producedSide.getAtom(j);
							if (eductAtom.getProperty(ATOM_ATOM_MAPPING) != null
									&& eductAtom.getProperty(ATOM_ATOM_MAPPING).equals(productAtom.getProperty(ATOM_ATOM_MAPPING))) {
								reaction.addMapping(
										builder.newInstance(IMapping.class, eductAtom, productAtom));
								mappingCount++;
								break;
							}
						}
					}
					logger.info("Mapped atom pairs: " + mappingCount);
					//						reactionSet.addReaction(reaction);
//				}
			}

		} catch (IllegalArgumentException | IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		if (reactionSet.getReactionCount() == 0) {
			reaction.setProperties(properties);
			reactionSet.addReaction(reaction);
			reactionCount++;
		}
		else if (reactionSet.getReaction(reactionSet.getReactionCount()-1).getID() != reaction.getID()) {
			reaction.setProperties(properties);
			reactionSet.addReaction(reaction);
			reactionCount++;
		}
		
		return reactionSet;
	}

	private double calculateConstant(String input) {
		double value = 0;
		if (input.contains("|")) {
			String[] arr1 = input.split("\\|");
			int cpt = 0;
			for (int i = 0; i < arr1.length; i++) {
				String val = arr1[i];
				if (val.contains(" - ")) {
					String[] arr2 = val.split(" - ");
					int cpt2 = 0;
					float val2 = 0;
					for (int j = 0; j < arr2.length; j++) {
						val2 += Double.parseDouble(arr2[j]);
						cpt2++;
					}
					val = String.valueOf(val2/cpt2);
				}
				value += Double.parseDouble(arr1[i]);
				cpt++;
			}
			return value/cpt;
		}
		else {
			String[] arr = input.split(" - ");
			int cpt = 0;
			for (int i = 0; i < arr.length; i++) {
				value += Double.parseDouble(arr[i]);
				cpt++;
			}
			return value/cpt;
		}
	}

	private void writeRDFChunk(IReactionSet reactions, String reactionFile) throws IOException {
		File directory = new File(new File(reactionFile).getParent(), "tempRDFiles/");
	    if (!directory.exists()){
	        directory.mkdir();
	    }
		File file = new File(directory.toString(), "chunk_" + files.size() + ".rdf");
		try (MDLV2000RDFWriter writer = new MDLV2000RDFWriter(new FileWriter(file))) {
			writer.write(reactions);
			//writer.setRdFieldsProperties(map);
			writer.close();
		} catch (CDKException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		files.add(file);
	}
	
	@Override
	public void close() throws IOException {
		input.close();
	}
	
	public boolean isAromatize() {
		return aromatize;
	}

	public void setAromatize(boolean aromatize) {
		this.aromatize = aromatize;
	}

	public boolean isSplit() {
		return split;
	}

	public void setSplit(boolean split) {
		this.split = split;
	}

	public String isFilename() {
		return filename;
	}

	public void setFilename(String filename) {
		this.filename = filename;
	}

	public List<File> getFiles() {
		return files;
	}

	public void setFiles(List<File> files) {
		this.files = files;
	}

	public int getReactionCount() {
		return reactionCount;
	}

	public void setReactionCount(int reactionCount) {
		this.reactionCount = reactionCount;
	}
	
}
