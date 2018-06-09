package ca.uwaterloo.eclipse.refactoring;

import java.io.File;
import java.io.IOException;
import java.lang.reflect.Field;
import java.lang.reflect.Method;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import ca.uwaterloo.eclipse.refactoring.coverage.CoverageStatus;
import ca.uwaterloo.eclipse.refactoring.coverage.LineCoverage;
import ca.uwaterloo.eclipse.refactoring.coverage.TestReport;
import ca.uwaterloo.eclipse.refactoring.coverage.TestReportResults;
import ca.uwaterloo.eclipse.refactoring.utility.IOHelper;
import ca.uwaterloo.eclipse.refactoring.utility.SourceDirectoryUtility;
import org.eclipse.core.resources.IFolder;
import org.eclipse.core.resources.IResource;
import org.eclipse.core.resources.IncrementalProjectBuilder;
import org.eclipse.core.runtime.CoreException;
import org.eclipse.core.runtime.IProgressMonitor;
import org.eclipse.core.runtime.NullProgressMonitor;
import org.eclipse.core.runtime.Path;
import org.eclipse.debug.core.DebugPlugin;
import org.eclipse.debug.core.ILaunch;
import org.eclipse.debug.core.ILaunchConfiguration;
import org.eclipse.debug.core.ILaunchConfigurationType;
import org.eclipse.debug.core.ILaunchConfigurationWorkingCopy;
import org.eclipse.debug.core.ILaunchManager;
import org.eclipse.debug.core.Launch;
import org.eclipse.jdt.core.IClasspathEntry;
import org.eclipse.jdt.core.ICompilationUnit;
import org.eclipse.jdt.core.IJavaProject;
import org.eclipse.jdt.core.IPackageFragment;
import org.eclipse.jdt.core.ISourceRange;
import org.eclipse.jdt.core.JavaCore;
import org.eclipse.jdt.core.JavaModelException;
import org.eclipse.jdt.core.ToolFactory;
import org.eclipse.jdt.core.formatter.CodeFormatter;
import org.eclipse.jdt.launching.IJavaLaunchConfigurationConstants;
import org.eclipse.jdt.launching.JavaLaunchDelegate;
import org.eclipse.text.edits.TextEdit;

public class ApplicationRunner extends JavaLaunchDelegate {
	public enum TestReportFileType {

		ORIGINAL(ORIGINAL_TEST_REPORT_CSV_FILE_NAME),
		AFTER_REFACTORING(TEST_REPORT_AFTER_REFACTORING_CSV);

		private String fileName;

		private TestReportFileType(String fileName) {
			this.fileName = fileName;
		}

		@Override
		public String toString() {
			return fileName;
		}
	}

	private static final String COVERAGE_REPORT_CSV_FILE_NAME = "coverage_report.csv";
	private static final String ORIGINAL_TEST_REPORT_CSV_FILE_NAME = "test_report.csv";
	private static final String TEST_REPORT_AFTER_REFACTORING_CSV = "test_report_offset.csv";
	private static final String JUNIT_LAUNCH_CONFIGURATION_TYPE = "org.eclipse.jdt.junit.launchconfig";
	private static final String JACOCO_EXEC_FILE = "/jacoco.exec";
	private static final String JACOCO_AGENT = "/lib/jacocoagent.jar";
	private static final String COVERAGE_JAR = "/lib/code-coverage-0.1-jar-with-dependencies.jar";
	private static final String LAUNCH_MAIN_TYPE_NAME_TESTS = "TestRunner";
	private static final String LAUNCH_MAIN_TYPE_NAME_COVERAGE_REPORT = "CoverageRunner";

	private final IJavaProject jProject;
	private final String classFolder;
	private final String reportFolder;
	private ILaunchConfiguration launchConfigurationForTest;
	private ILaunchConfiguration launchConfigurationForCoverage;
	private ILaunchManager launchManager;
	private Launch launchInstanceForTest;
	private Launch launchInstanceForCoverage;

	public ApplicationRunner(IJavaProject jProject, String classFolder, String reportFolder) throws IOException {
		launchManager = getLaunchManager();
		this.jProject = jProject;
		this.classFolder = classFolder;
		this.reportFolder = reportFolder;
	}

	public void launchApplication() throws CoreException, IOException {
		launchConfigurationForTest = null;
		launchConfigurationForCoverage = null;

		deleteLaunchConfiguration(launchManager);
		if (launchConfigurationForTest == null) {
			createBuildPath();
			createRunnableClassForTest();
			createRunnableClassForCoverage();
			createLaunchConfigurationForTest(launchManager);
			createLaunchConfigurationForCoverage(launchManager);
		}
		jProject.getProject().refreshLocal(IResource.DEPTH_INFINITE, new NullProgressMonitor());
		jProject.getProject().build(IncrementalProjectBuilder.INCREMENTAL_BUILD, new NullProgressMonitor());

		launchTest();
		launchCoverage();
	}

	private void cleanupLaunchManager(ILaunch launch){
		try {
			clearFieldByInvokingClear(launchManager.getClass().
			        getDeclaredField("fListeners"),launchManager);
			clearFieldByInvokingClear( DebugPlugin.getDefault().getClass().
			        getDeclaredField("fEventQueue"), DebugPlugin.getDefault());
			clearFieldByInvokingClear(launch.getClass().
			        getDeclaredField("fProcesses"),launch);
		} catch (NoSuchFieldException | SecurityException e) {
			e.printStackTrace();
		}
		
//		Field fListener =launchManager.getClass().
//		        getDeclaredField("fListeners");
//		fListener.setAccessible(true);
//		Object fListenerValue=fListener.get(launchManager);
//		Method clearListener=fListenerValue.getClass().getMethod("clear",new Class[]{});
//		clearListener.invoke(fListenerValue);
//		
//		Field fEventQueue = DebugPlugin.getDefault().getClass().
//		        getDeclaredField("fEventQueue");
//		fEventQueue.setAccessible(true);
//		Object fEventQueueValue=fEventQueue.get(DebugPlugin.getDefault());
//		Method clearEventQueue=fEventQueueValue.getClass().getMethod("clear",new Class[]{});
//		clearEventQueue.invoke(fEventQueueValue);
//		
//		Field fProcesses = launch.getClass().
//		        getDeclaredField("fProcesses");
//		fProcesses.setAccessible(true);
//		Object fProcessesValue=fProcesses.get(launch);
//		Method clearProcesses=fProcessesValue.getClass().getMethod("clear",new Class[]{});
//		clearProcesses.invoke(fProcessesValue);
	}

	private void clearFieldByInvokingClear(Field field, Object receiver) {
		try {
			field.setAccessible(true);
			Object fieldObject=field.get(receiver);
			Method method=fieldObject.getClass().getMethod("clear",new Class[]{});
			method.invoke(fieldObject);
		} catch (Exception e){
			e.printStackTrace();
		}
	}

	private void deleteLaunchConfiguration(ILaunchManager launchManager) throws CoreException {
		for (ILaunchConfiguration configuration : launchManager.getLaunchConfigurations()) {
			if (configuration.getName().equals(getRunConfigurationNameForTest())) {
				configuration.delete();
			}
			if (configuration.getName().equals(getRunConfigurationNameForCoverage())) {
				configuration.delete();
			}
		}
	}

	private ILaunchConfiguration getLaunchConfiguration(ILaunchManager launchManager, String launchConfiguratioName) throws CoreException {
		for (ILaunchConfiguration configuration : launchManager.getLaunchConfigurations()) {
			if (configuration.getName().equals(launchConfiguratioName)) {
				return configuration;
			}
		}
		return null;
	}

	public void launchCoverage() throws CoreException {
		if (launchConfigurationForCoverage==null)
			launchConfigurationForCoverage = getLaunchConfiguration(launchManager, getRunConfigurationNameForCoverage());
		if (launchInstanceForCoverage==null)
			launchInstanceForCoverage = new Launch(launchConfigurationForCoverage, ILaunchManager.RUN_MODE, null);
		launch(launchConfigurationForCoverage, ILaunchManager.RUN_MODE, launchInstanceForCoverage, new NullProgressMonitor());
		
		while (!launchInstanceForCoverage.isTerminated()) {
			try {
				Thread.sleep(100);
			} catch (InterruptedException e) {
				e.printStackTrace();
			}
		}
		//launchInstanceForCoverage.terminate();
		launchConfigurationForCoverage=null;
		launchManager.removeLaunch(launchInstanceForCoverage);
		cleanupLaunchManager(launchInstanceForCoverage);
	}

	public void launchTest() throws CoreException {
		if (launchConfigurationForTest==null)
			launchConfigurationForTest = getLaunchConfiguration(launchManager, getRunConfigurationNameForTest());
		if (launchInstanceForTest==null)
			launchInstanceForTest = new Launch(launchConfigurationForTest, ILaunchManager.RUN_MODE, null);
		launch(launchConfigurationForTest, ILaunchManager.RUN_MODE, launchInstanceForTest, new NullProgressMonitor());

		while (!launchInstanceForTest.isTerminated()) {
			try {
				Thread.sleep(100);
			} catch (InterruptedException e) {
				e.printStackTrace();
			}
		}
		launchConfigurationForTest=null;
		launchManager.removeLaunch(launchInstanceForTest);
		cleanupLaunchManager(launchInstanceForTest);
	}

	private void createBuildPath() throws JavaModelException, IOException {
		IClasspathEntry[] entries = jProject.getRawClasspath();
		String libraryPath = new Path(new File(ApplicationRunner.class.getProtectionDomain().getCodeSource().getLocation().getPath() + COVERAGE_JAR).toString()).toPortableString();

		String pathToSrcFolder = "/" + jProject.getProject().getName() + "/src";

		IClasspathEntry[] newEntries = null;
		if (!checkClassPathExist(entries, pathToSrcFolder)) {
			newEntries = recreatePathEntries(entries, pathToSrcFolder, JavaCore.newSourceEntry(new Path(pathToSrcFolder)));
		}
		if (!checkClassPathExist(entries, libraryPath)) {
			newEntries = newEntries == null ? entries : newEntries;
			newEntries = recreatePathEntries(newEntries, pathToSrcFolder, JavaCore.newLibraryEntry(new Path(libraryPath), null, null, false));
		}
		if (newEntries != null)
			jProject.setRawClasspath(newEntries, null);
	}

	private IClasspathEntry[] recreatePathEntries(IClasspathEntry[] entries, String pathToSrcFolder, IClasspathEntry entry) {
		IClasspathEntry[] newEntries;
		newEntries = new IClasspathEntry[entries.length + 1];
		System.arraycopy(entries, 0, newEntries, 0, entries.length);
		newEntries[newEntries.length - 1] = entry;
		return newEntries;
	}

	private boolean checkClassPathExist(IClasspathEntry[] entries, String libraryPath) {
		for (IClasspathEntry classEntry : entries) {
			if (classEntry.getPath().toPortableString().equals(libraryPath))
				return true;
		}
		return false;
	}

	private String getRunConfigurationNameForTest() {
		return jProject.getProject().getName() + "-test";
	}

	private String getRunConfigurationNameForCoverage() {
		return jProject.getProject().getName() + "-coverage";
	}

	/**
	 * Creates Application launch configuration to run unit tests
	 * 
	 * @param launchManager
	 * @return
	 * @throws CoreException
	 */
	private ILaunchConfiguration createLaunchConfigurationForTest(ILaunchManager launchManager) throws CoreException {
		Map<String, String> attributes = new HashMap<>();
		attributes.put(IJavaLaunchConfigurationConstants.ATTR_PROJECT_NAME, jProject.getProject().getName());
		attributes.put(IJavaLaunchConfigurationConstants.ATTR_MAIN_TYPE_NAME, LAUNCH_MAIN_TYPE_NAME_TESTS);
		attributes.put(IJavaLaunchConfigurationConstants.ATTR_PROGRAM_ARGUMENTS, getProgramArguments());

		attributes.put(IJavaLaunchConfigurationConstants.ATTR_VM_ARGUMENTS, getVMArguments());
		return createLaunchConfiguration(launchManager, getRunConfigurationNameForTest(), IJavaLaunchConfigurationConstants.ID_JAVA_APPLICATION, attributes);
	}

	/**
	 * Creates JUnit launch configuration.
	 * 
	 * @param launchManager
	 * @return
	 * @throws CoreException
	 */
	@SuppressWarnings("unused")
	private ILaunchConfiguration createLaunchConfigurationForJUnitTest(ILaunchManager launchManager) throws CoreException {
		Map<String, String> attributes = new HashMap<>();
		attributes.put("org.eclipse.jdt.launching.MAIN_TYPE", LAUNCH_MAIN_TYPE_NAME_TESTS);
		attributes.put("org.eclipse.jdt.launching.PROJECT_ATTR", jProject.getProject().getName());
		attributes.put(IJavaLaunchConfigurationConstants.ATTR_VM_ARGUMENTS, getVMArguments());
		attributes.put("org.eclipse.jdt.junit.TEST_KIND", "org.eclipse.jdt.junit.loader.junit4");

		return createLaunchConfiguration(launchManager, getRunConfigurationNameForTest(), JUNIT_LAUNCH_CONFIGURATION_TYPE, attributes);
	}

	private ILaunchConfiguration createLaunchConfigurationForCoverage(ILaunchManager launchManager) throws CoreException {
		Map<String, String> attributes = new HashMap<>();
		attributes.put(IJavaLaunchConfigurationConstants.ATTR_PROJECT_NAME, jProject.getProject().getName());
		attributes.put(IJavaLaunchConfigurationConstants.ATTR_MAIN_TYPE_NAME, LAUNCH_MAIN_TYPE_NAME_COVERAGE_REPORT);
		attributes.put(IJavaLaunchConfigurationConstants.ATTR_PROGRAM_ARGUMENTS, getProgramArguments());
		return createLaunchConfiguration(launchManager, getRunConfigurationNameForCoverage(), IJavaLaunchConfigurationConstants.ID_JAVA_APPLICATION, attributes);
	}

	private ILaunchConfiguration createLaunchConfiguration(ILaunchManager launchManager, String runConfigurationName, String launchConfigurationType, Map<String, String> attribtues) throws CoreException {
		ILaunchConfigurationType launchType = launchManager.getLaunchConfigurationType(launchConfigurationType);

		ILaunchConfigurationWorkingCopy workingCopy = launchType.newInstance(null, runConfigurationName);
		List<IResource> resources = new ArrayList<IResource>();
		resources.add(jProject.getProject());
		IResource[] resourcesArray = new IResource[resources.size()];
		resourcesArray = resources.toArray(resourcesArray);
		workingCopy.setMappedResources(resourcesArray);
		workingCopy.setAttributes(attribtues);
		return workingCopy.doSave();
	}

	private String getProgramArguments() throws JavaModelException {
		StringBuffer programArguments = new StringBuffer();
		programArguments.append("-p ");
		programArguments.append("\"" + jProject.getProject().getLocation() + "\"");
		programArguments.append(" ");
		programArguments.append("-s ");
		programArguments.append("\"" + String.join(",", SourceDirectoryUtility.getAllSourceDirectories(jProject)) + "\"");
		programArguments.append(" ");
		programArguments.append("-c ");
		programArguments.append("\"" + classFolder + "\"");
		programArguments.append(" ");
		programArguments.append("-r ");
		programArguments.append("\"" + reportFolder + "\"");
		return programArguments.toString();
	}

	private String getVMArguments() {
		StringBuffer vmArguments = new StringBuffer();
		try {
			vmArguments.append("-javaagent:");
			vmArguments.append(new File(ApplicationRunner.class.getProtectionDomain().getCodeSource().getLocation().getPath() + JACOCO_AGENT).getCanonicalFile().toString());
			vmArguments.append("=");
			vmArguments.append("destfile");
			vmArguments.append("=");
			vmArguments.append(jProject.getProject().getLocation());
			vmArguments.append(JACOCO_EXEC_FILE);
			// vmArguments.append("=");
			// vmArguments.append("-D" + REPORT_DIRECTORY_SYSTEM_PROPERTY_KEY);
			// vmArguments.append("=");
			// vmArguments.append("\"" + reportFolder + "\"");

		} catch (IOException e) {
			e.printStackTrace();
		}
		return vmArguments.toString();
	}

	private void createRunnableClassForTest() throws CoreException, JavaModelException {
		createRunnableClass(LAUNCH_MAIN_TYPE_NAME_TESTS, "TestAndCoverageEngine.run(args);");

		/*
		 * This code should be added to TestAndCoverageEngine when we create
		 * JUnit Run Configuration
		 * 
		 * @RunWith(TestRunner.AllTestsRunner.class) public final class
		 * TestRunner { public static class AllTestsRunner extends Suite {
		 * public AllTestsRunner(final Class<?> clazz) throws
		 * InitializationError { super( clazz, TestUtility .findClasses(
		 * "/Users/Shahriar/Documents/Workspace/jdeodorant-workspace/TestData/jruby/jruby-1.4.0/bin"
		 * )); System.out.println(System.getenv("binPath")); } } }
		 */
	}

	private void createRunnableClassForCoverage() throws CoreException, JavaModelException {
		createRunnableClass(LAUNCH_MAIN_TYPE_NAME_COVERAGE_REPORT, "TestAndCoverageEngine.generateReport(args);");
	}

	public void formatSourceCode(ICompilationUnit compilationUnit, IProgressMonitor monitor) throws JavaModelException {
		compilationUnit.becomeWorkingCopy(new NullProgressMonitor());
		CodeFormatter formatter = ToolFactory.createCodeFormatter(null);
		ISourceRange range = compilationUnit.getSourceRange();
		TextEdit formatEdit = formatter.format(CodeFormatter.K_COMPILATION_UNIT, compilationUnit.getSource(), range.getOffset(), range.getLength(), 0, null);
		if (formatEdit != null && formatEdit.hasChildren()) {
			compilationUnit.applyTextEdit(formatEdit, monitor);
		} else {
			monitor.done();
		}
		compilationUnit.commitWorkingCopy(true, new NullProgressMonitor());
	}

	private void createRunnableClass(String className, String statement) throws CoreException, JavaModelException {
		IFolder sourceFolder = jProject.getProject().getFolder("src");
		if (!sourceFolder.exists())
			sourceFolder.create(true, true, new NullProgressMonitor());
		IPackageFragment defaultPackage = jProject.getPackageFragmentRoot(sourceFolder).createPackageFragment("", true, null);
		StringBuffer buffer = new StringBuffer();
		buffer.append("import ca.concordia.jdeodorant.coverage.tools.TestAndCoverageEngine;");
		buffer.append("\n");
		buffer.append("public class " + className + "{");
		buffer.append("\n");
		buffer.append("public static void main(String[] args) {");
		buffer.append("\n");
		buffer.append(statement);
		buffer.append("\n");
		buffer.append("}");
		buffer.append("\n");
		buffer.append("}");
		String compilationUnitName = className + ".java";
		ICompilationUnit compilationUnit = defaultPackage.getCompilationUnit(compilationUnitName);
		if (compilationUnit.exists())
			compilationUnit.delete(true, new NullProgressMonitor());

		compilationUnit = defaultPackage.createCompilationUnit(compilationUnitName, buffer.toString(), true, new NullProgressMonitor());
		formatSourceCode(compilationUnit, new NullProgressMonitor());
	}

	public static List<LineCoverage> readCoverageFile(String pathToReportFolder) {
		List<LineCoverage> lineCoverageList = new ArrayList<>();
		List<String> lines = IOHelper.readFileByLine(pathToReportFolder + "/" + COVERAGE_REPORT_CSV_FILE_NAME);
		for (int i = 1; i < lines.size(); i++) {
			String line = lines.get(i);
			String[] values = line.split(",");
			lineCoverageList.add(new LineCoverage(values[0], values[1], Integer.parseInt(values[2]), CoverageStatus.valueOf(values[3])));
		}
		return lineCoverageList;
	}

	public static TestReportResults readTestFile(String pathToReportFolder, TestReportFileType fileType) {
		TestReportResults testReportList = new TestReportResults();
		List<String> lines = IOHelper.readFileByLine(pathToReportFolder + "/" + fileType);
		for (int i = 1; i < lines.size(); i++) {
			String line = lines.get(i);
			String[] values = line.split(",");
			testReportList.addTestResult(new TestReport(values[0], values[1], values[2]));
		}
		return testReportList;
	}
}
