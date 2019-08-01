package ca.uwaterloo.jrefactoring.match;

import org.eclipse.core.filebuffers.FileBuffers;
import org.eclipse.core.filebuffers.ITextFileBuffer;
import org.eclipse.core.filebuffers.ITextFileBufferManager;
import org.eclipse.core.filebuffers.LocationKind;
import org.eclipse.core.runtime.CoreException;
import org.eclipse.core.runtime.IPath;
import org.eclipse.jdt.core.*;
import org.eclipse.jface.text.BadLocationException;
import org.eclipse.jface.text.IDocument;

import ca.concordia.jdeodorant.clone.parsers.CloneGroup;

import java.util.*;

public class CloneItem {
	private CloneGroupInfo cloneGroupBelonging;
	private int numberOfPDGNodes = 0;
	private int numberOfCloneStatements = 0;
	private int startOffsetOfCodeFragment = 0;
	private int endOffsetOfCodeFragment = 0;
	private String sourceCode = null;
	private IMethod iMethod;
	private int cloneFragmentID;
	private ICompilationUnit iCompilationUnit;
	private String sourceFolder;
	private String cloneClass, clonePackage;
	
	public CloneGroupInfo getCloneGroupBelonging() {
		return cloneGroupBelonging;
	}

	public void setCloneGroupBelonging(CloneGroupInfo cloneGroupBelonging) {
		this.cloneGroupBelonging = cloneGroupBelonging;
	}
	
	public int getNumberOfPDGNodes() {
		return numberOfPDGNodes;
	}

	public void setNumberOfPDGNodes(int numberOfPDGNodes) {
		this.numberOfPDGNodes = numberOfPDGNodes;
	}
	
	public int getNumberOfCloneStatements() {
		return numberOfCloneStatements;
	}

	public void setNumberOfCloneStatements(int numberOfCloneStatements) {
		this.numberOfCloneStatements = numberOfCloneStatements;
	}
	
	public int getStartOffsetOfCodeFragment() {
		return startOffsetOfCodeFragment;
	}

	public void setStartOffsetOfCodeFragment(int startOffsetOfCodeFragment) {
		this.startOffsetOfCodeFragment = startOffsetOfCodeFragment;
		this.sourceCode = null;
	}
	
	public int getEndOffsetOfCodeFragment() {
		return endOffsetOfCodeFragment;
	}

	public void setEndOffsetOfCodeFragment(int endOffsetOfCodeFragment) {
		this.endOffsetOfCodeFragment = endOffsetOfCodeFragment;
	}
	
	public String getSourceCode() {
		if (sourceCode == null) {
			try {
				sourceCode = "";
				if (startOffsetOfCodeFragment < Integer.MAX_VALUE && endOffsetOfCodeFragment > -1)
					sourceCode = getSourceCodeStringFromICompilationUnit(startOffsetOfCodeFragment, endOffsetOfCodeFragment, iCompilationUnit);
			} catch (BadLocationException e) {
				e.printStackTrace();
			}
		}
		return sourceCode;
	}
	
	public String getFirstMethodSignature() {
		if (iMethod != null)
			return getMethodJavaSignature(iMethod);
		return "";
	}

	public void setContaingingIMethodFirst(IMethod iMethod) {
		this.iMethod = iMethod;
		this.iCompilationUnit = iMethod.getCompilationUnit();
	}
	
	public int getCloneFragmentID() {
		return cloneFragmentID;
	}

	public void setCloneFragmentID(int cloneFragmentID) {
		this.cloneFragmentID = cloneFragmentID;
	}
	
	private String getSourceCodeStringFromICompilationUnit(int startOffset, int endOffset, ICompilationUnit iCompilationUnit) throws BadLocationException {
		ITextFileBufferManager bufferManager = FileBuffers.getTextFileBufferManager();
		IPath iPath = iCompilationUnit.getPath();

		try {
			bufferManager.connect(iPath, LocationKind.IFILE, null);
			ITextFileBuffer textFileBuffer = bufferManager.getTextFileBuffer(iPath, LocationKind.IFILE);
			IDocument iDocument = textFileBuffer.getDocument();
			return iDocument.get(startOffset, endOffset - startOffset);

		} catch (CoreException e) {
			e.printStackTrace();
		}
		return null;
	}
	
	public String getContainingFileFirst() {
		return iCompilationUnit.getPath().toPortableString();
	}
	
	public ICompilationUnit getICompilationUnit() {
		return iCompilationUnit;
	}

	public void setICompilationUnitFirst(ICompilationUnit iCompilationUnit) {
		this.iCompilationUnit = iCompilationUnit;
	}
	
	public void setSourceFolder(String sourceFolder) {
		this.sourceFolder = sourceFolder;
	}

	public String getSourceFolder() {
		return this.sourceFolder;
	}
	
	public String getCloneClass() {
		return this.cloneClass;
	}

	public void setCloneClass(String cloneClass) {
		this.cloneClass = cloneClass;
	}
	
	public String getClonePackage() {
		return this.clonePackage;
	}

	public void setClonePackage(String clonePackage) {
		this.clonePackage = clonePackage;
	}
	
	private String getMethodJavaSignature(IMethod iMethod) {

		StringBuilder toReturn = new StringBuilder();

		//toReturn.append(iMethod.getDeclaringType().getFullyQualifiedName());
		try {
			toReturn.append(Signature.toString(iMethod.getReturnType()));
		} catch (IllegalArgumentException  | JavaModelException e) {

		}
		toReturn.append(" ");
		toReturn.append(iMethod.getElementName());
		toReturn.append("(");

		String comma = "";
		for (String type : iMethod .getParameterTypes()) {
			toReturn.append(comma);
			comma = ", ";
			toReturn.append(Signature.toString(type));
		}
		toReturn.append(")");

		return toReturn.toString();
	}
	
}
