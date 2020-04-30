package com.nih.tools;

import org.openscience.cdk.AtomContainer;
import org.openscience.cdk.interfaces.IAtomContainer;

public class Clone {
	
	public static synchronized IAtomContainer deepClone(IAtomContainer ac) throws CloneNotSupportedException {
            IAtomContainer acClone = new AtomContainer(ac).clone();
            /*Set IDs as CDK clone doesn't*/
            for (int i = 0; i < ac.getAtomCount(); i++) {
                acClone.getAtom(i).setID(ac.getAtom(i).getID());
            }
            for (int i = 0; i < ac.getBondCount(); i++) {
                acClone.getBond(i).setID(ac.getBond(i).getID());
            }
            acClone.setID(ac.getID());
            acClone.addProperties(ac.getProperties());

        return acClone;
    }
}
