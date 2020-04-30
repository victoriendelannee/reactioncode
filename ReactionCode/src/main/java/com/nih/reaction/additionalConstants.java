package com.nih.reaction;

public class additionalConstants {

	public final static String ATOMCONTAINER_INDEX = "IAtomContainerIndex";

	public final static String BOND_BOND_MAPPING = "BondBondMapping";
	
	public final static String REACTING_CENTER_STATUS = "ReactingCenterStatus";

	/** Flag is set if bond is a stereocenter. */
	public final static int IS_STEREOCENTER = 0x0001;

	/** Flag is set if atom or bond is part of the reation center. */
	public final static int REACTION_CENTER = 0x0002;

	/** Flag is set if atom is not in product. */
	public final static int LEAVING_ATOM = 0x0003;

	/** Flag is set if bond is formed. */
	public final static int BOND_FORMED = 0x0004;

	/** Flag is set if bond is cleaved. */
	public final static int BOND_CLEAVED = 0x0005;

	/** Flag is set if bond has an order change. */
	public final static int BOND_ORDER = 0x0006;

	/** Flag is set if bond has a stereo Change. */
	public final static int BOND_STEREO = 0x0007;
	
	/** Flag is set if an atom gain or loose a charge. */
	public final static int CHARGE_CHANGE = 0x0008;
	
	/** get the type of bond cleaved or formed */
	public final static String BOND_CHANGE_INFORMATION = "bondChangeInformation";
	
	/** get the order Change form previous state to new state */
	public final static String BOND_ORDER_CHANGE = "BondOrderChange";

	/** get the stereo Change form previous state to new state */
	public final static String BOND_STEREO_CHANGE = "BondStereoChange";
	
	/** type of stereo (tetra, planar,...)*/
	public final static String STEREO_TYPE = "StereoType";
	
	/** symbol of stereo (DB1, OH1, ...*/
	public final static String STEREO_SYMBOL = "StereoSymbol";
	
	/** symbol of stereo (DB1, OH1, ...*/
	public final static String STEREO_SHORTHAND = "StereoShorthand";
	
	/** get the stereo Change form previous state to new state */
	public final static String ATOM_CHARGE_CHANGE = "AtomChargeChange";
	
	public final static int BOND_ORDER_GAIN = 0x0009;
	
	public final static int BOND_ORDER_REDUCED = 0x000A;
	
//	public final static int SINGLE_TO_DOUBLE = 0x0007;
//	
//	public final static int SINGLE_TO_TRIPLE = 0x0104;
//	
//	public final static int SINGLE_TO_QUADRIPLE = 0x0100;
//	
//	public final static int DOUBLE_TO_SINGLE = 0x0008;
//	
//	public final static int DOUBLE_TO_TRIPLE = 0x0009;
//	
//	public final static int DOUBLE_TO_QUDRIPLE= 0x0101;
//	
//	public final static int TRIPLE_TO_DOUBLE = 0x0010;
//	
//	public final static int TRIPLE_TO_QUADRIPLE = 0x0102;
//	
//	public final static int QUADRIPLE_TO_TRIPLE = 0x0103;
	
	public final static int DOWN_TO_UP = 0x000B;
	
	public final static int UP_TO_DOWN = 0x000C;
	
	public final static int DOWN_INVERTED_TO_UP_INVERTED = 0x000D;
	
	public final static int UP_INVERTED_TO_DOWN_INVERTED = 0x000E;
	
	public final static int E_TO_Z = 0x000F;
	
	public final static int Z_TO_E = 0x0010;
	
	public final static int NONE_TO_UP = 0x0011;
	
	public final static int NONE_TO_DOWN = 0x0012;
	
	public final static int NONE_TO_UP_INVERTED = 0x0013;
	
	public final static int NONE_TO_DOWN_INVERTED = 0x0014;
	
	public final static int NONE_TO_Z = 0x0015;
	
	public final static int NONE_TO_E = 0x0016;
	
	public final static int UP_TO_NONE = 0x0017;
	
	public final static int DOWN_TO_NONE = 0x0018;
	
	public final static int UP_INVERTED_TO_NONE = 0x0019;
	
	public final static int DOWN_INVERTED_TO_NONE = 0x001A;
	
	public final static int Z_TO_NONE = 0x001B;
	
	public final static int E_TO_NONE = 0x001C;
	
	public final static int DOWN__TO_UP_INVERTED = 0x001D;
	
	public final static int UP_TO_DOWN_INVERTED = 0x001E;
	
	public final static int DOWN_INVERTED_TO_UP = 0x001F;
	
	public final static int UP_INVERTED_TO_DOWN = 0x0020;
	
	public final static int GAIN_ONE_CHARGE = 0x0021;
	
	public final static int LOOSE_ONE_CHARGE = 0x0022;
	
}
