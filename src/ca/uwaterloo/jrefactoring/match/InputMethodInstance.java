package ca.uwaterloo.jrefactoring.match;

import gr.uom.java.ast.decomposition.cfg.PDG;
import org.eclipse.jdt.core.IMethod;

public class InputMethodInstance {
	private final IMethod iMethod;
	private final int startOffset;
	private final int endOffset;
	private final PDG pdg;
	
	public InputMethodInstance(IMethod iMethod, int startOffset, int endOffset, PDG pdg) {
		this.iMethod = iMethod;
		this.startOffset = startOffset;
		this.endOffset = endOffset;
		this.pdg = pdg;
	}
	
	public IMethod getIMethod() {
		return iMethod;
	}

	public int getStartOffset() {
		return startOffset;
	}

	public int getEndOffset() {
		return endOffset;
	}

	public PDG getPDG() {
		return pdg;
	}

}
