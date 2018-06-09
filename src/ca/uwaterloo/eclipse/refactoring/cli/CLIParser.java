package ca.uwaterloo.eclipse.refactoring.cli;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import ca.uwaterloo.eclipse.refactoring.parsers.CloneToolParserType;
import ca.uwaterloo.eclipse.refactoring.utility.Mailer;
import org.apache.commons.cli.BasicParser;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;


/**
 * 
 * @author Davood Mazinanian Creates the command line parser for this
 *         application, this is also a facade
 * 
 *         This otherargs option must include: <br />
 *         <ul>
 *         <li>When using CCFinder result files, this argument must include the
 *         path of the preprocessor files (path to the .ccfxprepdir folder)
 *         followed by the path to the src folder (relative to the project
 *         path), separated by comma. For exmaple, if project path is
 *         "path/to/project" and src folder is "path/to/project/src", one must
 *         provide "src".</li>
 *         <li>When using Nicad, this argument only includes the path to the src
 *         folder of the project, when used Nicad (because Nicad makes absolute
 *         paths)
 *         <li>
 *         </ul>
 */
public class CLIParser {

	private static Logger LOGGER = LoggerFactory.getLogger(CLIParser.class);

	public enum ApplicationMode {
		/** Only parse the output files of a clone tool */
		PARSE,
		/** Parse the output files of a clone tool, then start analysis */
		PARSE_AND_ANALYZE,
		/**
		 * Only analyze the existing excel file create from a parse execution of
		 * this tool using PARSE mode
		 */
		ANALYZE_EXISTING, Refactor,
	}

	private CommandLine cmdLine;
	private Options options;

	@SuppressWarnings("static-access")
	public CLIParser(String[] args) {

		options = new Options();

		options.addOption(OptionBuilder
				.withArgName("")
				.withLongOpt("mode")
				.withDescription("Mode of operation, one of: " + Arrays.toString(ApplicationMode.values()))
				.hasArg()
				.create("m"));

		options.addOption(OptionBuilder
				.withArgName("")
				.withLongOpt("project")
				.withDescription("Exact name of the project opened in the workspace")
				.hasArg()
				.create("p"));
		
		options.addOption(OptionBuilder
				.withArgName("")
				.withLongOpt("project-description")
				.withDescription("Path to the .project file of the project")
				.hasArg()
				.create("pd"));

		options.addOption(OptionBuilder
				.withArgName("")
				.withLongOpt("excelfile")
				.withDescription("Path to the excel file")
				.hasArg()
				.create("x"));

		options.addOption(OptionBuilder
				.withArgName("")
				.withLongOpt("tooloutputfile")
				.withDescription("Path to the output file of the clone tool")
				.hasArg()
				.create("i"));

		options.addOption(OptionBuilder
				.withArgName("")
				.withLongOpt("extra-args")
				.withDescription("Comma separated list of extra arguments. "
						+ "See documentation for more information.")
				.hasArgs()
				.withValueSeparator(',')
				.create("xargs"));

		options.addOption(OptionBuilder
				.withArgName("")
				.withLongOpt("tool")
				.withDescription("Name of the tool, one of: " + 
						Arrays.toString(CloneToolParserType.values()))
				.hasArg()
				.create("t"));

		options.addOption(OptionBuilder
				.withArgName("")
				.withLongOpt("help")
				.withDescription("Display this help")
				.hasOptionalArg()
				.create("?"));

		options.addOption(OptionBuilder
				.withArgName("")
				.withLongOpt("row-start-from")
				.withDescription("The row in the excel file from which the analysis should start from "
						+ "(must be the first row in the clone group)")
				.hasArg()
				.create("r"));

		options.addOption(OptionBuilder
				.withArgName("")
				.withLongOpt("append-results")
				.withDescription("Determine weather the result files (e.g. CSV files) "
						+ "must be appended to the existing files. By default, the existing "
						+ "files will be overridden.")
				.hasOptionalArg()
				.create("a"));

		options.addOption(OptionBuilder
				.withArgName("")
				.withLongOpt("skip-groups")
				.withDescription("A comma separated list of clone groups in the excel file to be skipped.")
				.withValueSeparator(',')
				.hasArgs()
				.create("s"));

		options.addOption(OptionBuilder
				.withArgName("")
				.withLongOpt("test-packages")
				.withDescription("A comma-separated list of the name of test packages")
				.withValueSeparator(',')
				.hasArgs()
				.create("testpkgs"));

		options.addOption(OptionBuilder
				.withArgName("")
				.withLongOpt("test-source-folders")
				.withDescription("A comma-separated list of the name of test source folders")
				.withValueSeparator(',')
				.hasArgs()
				.create("testsrcs"));

		options.addOption(OptionBuilder
				.withArgName("")
				.withLongOpt("class-folder")
				.withDescription("Name of class folder which .class files are located")
				.hasArg()
				.create("cf"));
		
		options.addOption(OptionBuilder
				.withArgName("")
				.withLongOpt("run-tests")
				.withDescription("Run tests after applying each refactoring")
				.hasOptionalArg()
				.create("rt"));
		
		options.addOption(OptionBuilder
				.withArgName("")
				.withLongOpt("log-to-file")
				.withDescription("Create a log file from console output")
				.hasOptionalArg()
				.create("l"));
		
		options.addOption(OptionBuilder
				.withArgName("")
				.withLongOpt("group-ids")
				.withDescription("A comma-separated list of clone group IDs to be analyzed. "
						+ "Other clone groups in the file will be skipped.")
				.withValueSeparator(',')
				.hasArgs()
				.create("g"));

		options.addOption(OptionBuilder
				.withArgName("")
				.withLongOpt("debugging-enabled")
				.withDescription("Prevent Eclipse command-line tool to cancel jobs"
						+ "queued in Eclipse JobManager such as workbench job and etc.")
				.hasOptionalArgs()
				.create("de"));

		options.addOption(OptionBuilder
				.withArgName("")
				.withLongOpt("mail-server-ip")
				.withDescription("SMTP server IP address")
				.hasArg()
				.create("msrvr"));
		
		options.addOption(OptionBuilder
				.withArgName("")
				.withLongOpt("mail-server-port")
				.withDescription("SMTP server port number")
				.hasArg()
				.create("mport"));
		
		options.addOption(OptionBuilder
				.withArgName("")
				.withLongOpt("mail-server-security-type")
				.withDescription("Security type for mail server, one of the values: " + Arrays.toString(Mailer.SecurityType.values()))
				.hasArg()
				.create("msectype"));
		
		options.addOption(OptionBuilder
				.withArgName("")
				.withLongOpt("mail-server-authenticated")
				.withDescription("Is SMTP server authenticated")
				.hasOptionalArg()
				.create("mauth"));
		
		options.addOption(OptionBuilder
				.withArgName("")
				.withLongOpt("mail-server-user-name")
				.withDescription("SMTP user name")
				.hasArg()
				.create("muser"));
		
		options.addOption(OptionBuilder
				.withArgName("")
				.withLongOpt("mail-server-password")
				.withDescription("SMTP password")
				.hasArg()
				.create("mpass"));
		
		options.addOption(OptionBuilder
				.withArgName("")
				.withLongOpt("email-addresses")
				.withDescription("Address to which the analysis notifications should be sent")
				.withValueSeparator(',')
				.hasOptionalArgs()
				.create("em"));
		
		// create the Apache Commons CLI parser
		CommandLineParser parser = new BasicParser();
		try {
			// parse the command line arguments
			this.cmdLine = parser.parse(options, args);
		} catch (ParseException exp) {
			throw new RuntimeException("CLI parsing failed. Reason: "
					+ exp.getMessage());
		}
	}

	private String getValue(String key) {
		return getValue(key, null);
	}

	private String getValue(String key, String defaultValue) {
		return cmdLine.getOptionValue(key, defaultValue);
	}

	private boolean hasOption(String option) {
		return cmdLine.hasOption(option);
	}

	/**
	 * Get an array containing the ID of clone groups to be skipped from
	 * analysis
	 * 
	 * @return
	 */
	public int[] getCloneGroupIDsToSkip() {
		return getIntArrayFromStringArrayForOption("skip-groups");
	}

	private int[] getIntArrayFromStringArrayForOption(String option) {
		int[] toReturn = new int[] {};
		String[] skippedCloneGroupsStringArray = cmdLine
				.getOptionValues(option);
		if (skippedCloneGroupsStringArray != null) {
			toReturn = new int[skippedCloneGroupsStringArray.length];
			for (int i = 0; i < skippedCloneGroupsStringArray.length; i++) {
				if (!"".equals(skippedCloneGroupsStringArray[i].trim()))
					try {
						toReturn[i] = Integer
								.parseInt(skippedCloneGroupsStringArray[i]);
					} catch (NumberFormatException nfex) {
						LOGGER.warn("Invalid clone group ID "
								+ skippedCloneGroupsStringArray[i]);
					}
			}
		}
		return toReturn;
	}

	/**
	 * Get an array containing the names of the test packages specified by user
	 * 
	 * @return
	 */
	public String[] getTestPackages() {
		return cmdLine.getOptionValues("test-packages");
	}

	/**
	 * Get an array containing the names of the test source folders specified by
	 * user
	 * 
	 * @return
	 */
	public String[] getTestSourceFolders() {
		return cmdLine.getOptionValues("test-source-folders");
	}

	/**
	 * Returns the values for "otherargs" option as an array of Strings
	 * 
	 * @param key
	 * @return
	 */
	public String[] getOtherArgs() {
		List<String> toReturn = new ArrayList<>();
		if (cmdLine.getOptionValue("xargs") != null) {
			String[] otherArgs = cmdLine.getOptionValues("xargs");
			for (int i = 0; i < otherArgs.length; i++) {
				if (!"".equals(otherArgs[i])) {
					// Remove starting and ending quotation marks
					if (otherArgs[i].startsWith("\""))
						otherArgs[i] = otherArgs[i].substring(1);
					if (otherArgs[i].endsWith("\""))
						otherArgs[i] = otherArgs[i].substring(0,
								otherArgs[i].length() - 1);
					toReturn.add(otherArgs[i]);
				}
			}
		}

		return toReturn.toArray(new String[] {});
	}

	public void printHelp() {
		HelpFormatter formatter = new HelpFormatter();
		formatter.printHelp("ca.concordia.jdeodorant.eclipse.commandline", options);
	}

	public ApplicationMode getApplicationMode(ApplicationMode defaultMode) {
		return ApplicationMode.valueOf(getValue("mode", defaultMode.toString()).toUpperCase());
	}

	public String getProjectName() {
		return getValue("project");
	}

	public boolean showHelp() {
		return hasOption("help");
	}

	public String getExcelFilePath() {
		return getValue("excelfile");
	}

	public int getStartingRow() {
		String rowStartFromString = getValue("row-start-from");
		int startFrom = 2;
		if (rowStartFromString != null && !"".equals(rowStartFromString)) {
			try {
				startFrom = Integer.parseInt(rowStartFromString);
			} catch (NumberFormatException nfe) {
				LOGGER.warn("Specified starting row is not valid, 2 selected as default.");
			}
		}
		return startFrom;
	}

	public boolean getAppendResults() {
		return hasOption("a");
	}

	public String getCloneToolOutputFilePath() {
		return getValue("tooloutputfile");
	}

	public String getCloneToolName() {
		return getValue("tool");
	}

	public Boolean hasCoverageReport() {
		// hasOption checks both longOpt name and short-one. cr stands for
		// coverage-report
		if (hasOption("cr"))
			return true;
		return false;
	}

	public String getClassFolder() {
		return getValue("cf");
	}

	public boolean runTests() {
		return hasOption("rt");
	}

	public boolean hasLogToFile() {
		return hasOption("l");
	}
	
	public boolean isDebuggingEnabled() {
		return hasOption("de");
	}
	
	public int[] getCloneGroupIDsToAnalyze() {
		return getIntArrayFromStringArrayForOption("group-ids");
	}

	public int getSMTPServerPort() {
		try {
			return Integer.parseInt(getValue("mport", "25"));
		} catch (NumberFormatException nfex) {
			LOGGER.warn("Invalid port number, 25 is selected by default");
			return 25;
		}
	}
	
	public String getSMTPServerAddress() {
		return getValue("msrvr", "127.0.0.1");
	}

	public String[] getNotificationEmailAddresses() {
		String[] optionValues = cmdLine.getOptionValues("em");
		if (optionValues != null) {
			return optionValues;
		}
		return new String[] {};
	}

	public boolean isMailServerAuthenticated() {
		return hasOption("mauth");
	}

	public String getMailServerUserName() {
		return getValue("muser", "");
	}

	public String getMailServerPassword() {
		return getValue("mpass", "");
	}

	public Mailer.SecurityType getMailServerSecurtyType() {
		try {
			return Mailer.SecurityType.valueOf(getValue("msectype", "none").toUpperCase());
		} catch (IllegalArgumentException ex) {
			return Mailer.SecurityType.NONE;
		}
	}

	public String getProjectDescritionFile() {
		return getValue("project-description");
	}
}
