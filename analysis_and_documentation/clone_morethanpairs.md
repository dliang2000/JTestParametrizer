# Clone Documentation - More than pairs

## Clones of 5
### Lucene-6.6.6  (minT = 30)
// 1236	src/test	org.apache.lucene.codecs.lucene54	TestLucene54DocValuesFormat	testSortedSetVariableLengthBigVsStoredFields	()V	89	95
// 1236	src/test	org.apache.lucene.codecs.lucene54	TestLucene54DocValuesFormat	testSortedVariableLengthManyVsStoredFields	()V	113	119
// 1236	src/test	org.apache.lucene.codecs.lucene54	TestLucene54DocValuesFormat	testTermsEnumFixedWidth	()V	121	127
// 1236	src/test	org.apache.lucene.codecs.lucene54	TestLucene54DocValuesFormat	testTermsEnumVariableWidth	()V	129	135
// 1236	src/test	org.apache.lucene.codecs.lucene54	TestLucene54DocValuesFormat	testTermsEnumRandomMany	()V	137	143

```java
@Slow
public void testSortedSetVariableLengthBigVsStoredFields() throws Exception {
  int numIterations = atLeast(1);
  for (int i = 0; i < numIterations; i++) {
    doTestSortedSetVsStoredFields(atLeast(300), 1, 32766, 16, 100);
  }
}


@Nightly
public void testSortedVariableLengthManyVsStoredFields() throws Exception {
  int numIterations = atLeast(1);
  for (int i = 0; i < numIterations; i++) {
    doTestSortedVsStoredFields(TestUtil.nextInt(random(), 1024, 2049), 1, 500);
  }
}

@Slow
public void testTermsEnumFixedWidth() throws Exception {
  int numIterations = atLeast(1);
  for (int i = 0; i < numIterations; i++) {
    doTestTermsEnumRandom(TestUtil.nextInt(random(), 1025, 5121), 10, 10);
  }
}

@Slow
public void testTermsEnumVariableWidth() throws Exception {
  int numIterations = atLeast(1);
  for (int i = 0; i < numIterations; i++) {
    doTestTermsEnumRandom(TestUtil.nextInt(random(), 1025, 5121), 1, 500);
  }
}

@Nightly
public void testTermsEnumRandomMany() throws Exception {
  int numIterations = atLeast(1);
  for (int i = 0; i < numIterations; i++) {
    doTestTermsEnumRandom(TestUtil.nextInt(random(), 1025, 8121), 1, 500);
  }
}
```

### Solr-6.6.6 (minT = 30)

// The following one is also minT = 100
// 742	src/test	org.apache.solr.highlight	HighlighterTest	testMultiValueAnalysisHighlight	()V	303	325
// 742	src/test	org.apache.solr.highlight	HighlighterTest	testMultiValueBestFragmentHighlight	()V	327	346
// 742	src/test	org.apache.solr.highlight	HighlighterTest	testDefaultFieldHighlight	()V	364	384
// 742	src/test	org.apache.solr.highlight	HighlighterTest	testTwoFieldHighlight	()V	405	425
// 742	src/test	org.apache.solr.highlight	HighlighterTest	testTermVecHighlight	()V	159	180


```java
@Test
public void testMultiValueAnalysisHighlight() {

  // do summarization using re-analysis of the field
  HashMap<String,String> args = new HashMap<>();
  args.put("hl", "true");
  args.put("hl.fl", "textgap");
  args.put("df", "textgap");
  TestHarness.LocalRequestFactory sumLRF = h.getRequestFactory(
    "standard", 0, 200, args);

  assertU(adoc("textgap", "first entry hasnt queryword",
      "textgap", "second entry has queryword long",
      "id", "1"));
  assertU(commit());
  assertU(optimize());
  assertQ("Basic summarization",
          sumLRF.makeRequest("long"),
          "//lst[@name='highlighting']/lst[@name='1']",
          "//lst[@name='1']/arr[@name='textgap']/str"
          );

}

@Test
public void testMultiValueBestFragmentHighlight() {
  HashMap<String,String> args = new HashMap<>();
  args.put("hl", "true");
  args.put("hl.fl", "textgap");
  args.put("df", "textgap");
  TestHarness.LocalRequestFactory sumLRF = h.getRequestFactory(
      "standard", 0, 200, args);

  assertU(adoc("textgap", "first entry has one word foo",
      "textgap", "second entry has both words foo bar",
      "id", "1"));
  assertU(commit());
  assertU(optimize());
  assertQ("Best fragment summarization",
      sumLRF.makeRequest("foo bar"),
      "//lst[@name='highlighting']/lst[@name='1']",
      "//lst[@name='1']/arr[@name='textgap']/str[.=\'second entry has both words <em>foo</em> <em>bar</em>\']"
  );
}

@Test
public void testDefaultFieldHighlight() {

  // do summarization using re-analysis of the field
  HashMap<String,String> args = new HashMap<>();
  args.put("hl", "true");
  args.put("df", "t_text");
  args.put("hl.fl", "");
  TestHarness.LocalRequestFactory sumLRF = h.getRequestFactory(
    "standard", 0, 200, args);

  assertU(adoc("t_text", "a long day's night", "id", "1"));
  assertU(commit());
  assertU(optimize());
  assertQ("Basic summarization",
      sumLRF.makeRequest("long"),
      "//lst[@name='highlighting']/lst[@name='1']",
      "//lst[@name='1']/arr[@name='t_text']/str"
  );

}


@Test
public void testTwoFieldHighlight() {

  // do summarization using re-analysis of the field
  HashMap<String,String> args = new HashMap<>();
  args.put("hl", "true");
  args.put("hl.fl", "t_text tv_text");
  TestHarness.LocalRequestFactory sumLRF = h.getRequestFactory(
    "standard", 0, 200, args);

  assertU(adoc("t_text", "a long day's night", "id", "1",
               "tv_text", "a long night's day"));
  assertU(commit());
  assertU(optimize());
  assertQ("Basic summarization",
          sumLRF.makeRequest("t_text:long"),
          "//lst[@name='highlighting']/lst[@name='1']",
          "//lst[@name='1']/arr[@name='t_text']/str",
          "//lst[@name='1']/arr[@name='tv_text']/str"
          );
}

@Test
public void testTermVecHighlight() {

  // do summarization using term vectors
  HashMap<String,String> args = new HashMap<>();
  args.put("hl", "true");
  args.put("hl.fl", "tv_text");
  args.put("hl.snippets", "2");
  TestHarness.LocalRequestFactory sumLRF = h.getRequestFactory(
    "standard",0,200,args);

  assertU(adoc("tv_text", LONG_TEXT,
               "id", "1"));
  assertU(commit());
  assertU(optimize());
  assertQ("Basic summarization",
          sumLRF.makeRequest("tv_text:long"),
          "//lst[@name='highlighting']/lst[@name='1']",
          "//lst[@name='1']/arr[@name='tv_text']/str[.='a <em>long</em> days night this should be a piece of text which']",
          "//arr[@name='tv_text']/str[.=' <em>long</em> fragments.']"
          );
}
```

// 35	src/test	org.apache.solr.handler.component	TestHttpShardHandlerFactory	testWhitelistHostCheckerSingleHost	()V	124	128
// 35	src/test	org.apache.solr.handler.component	TestHttpShardHandlerFactory	testWhitelistHostCheckerMultipleHost	()V	130	134
// 35	src/test	org.apache.solr.handler.component	TestHttpShardHandlerFactory	testWhitelistHostCheckerNoProtocolInParameter	()V	142	146
// 35	src/test	org.apache.solr.handler.component	TestHttpShardHandlerFactory	testWhitelistHostCheckerNonWhitelistedHostHttps	()V	172	176
// 35	src/test	org.apache.solr.handler.component	TestHttpShardHandlerFactory	testWhitelistHostCheckerCoreSpecific	()V	190	195

```java
@Test
public void testWhitelistHostCheckerSingleHost() {
  WhitelistHostChecker checker = new WhitelistHostChecker("http://abc-1.com:8983/solr", true);
  checker.checkWhitelist("http://abc-1.com:8983/solr", Arrays.asList(new String[]{"http://abc-1.com:8983/solr"}));
}

@Test
public void testWhitelistHostCheckerMultipleHost() {
  WhitelistHostChecker checker = new WhitelistHostChecker("http://abc-1.com:8983, http://abc-2.com:8983, http://abc-3.com:8983", true);
  checker.checkWhitelist("http://abc-1.com:8983/solr", Arrays.asList(new String[]{"http://abc-1.com:8983/solr"}));
}

@Test
public void testWhitelistHostCheckerNoProtocolInParameter() {
  WhitelistHostChecker checker = new WhitelistHostChecker("http://abc-1.com:8983, http://abc-2.com:8983, http://abc-3.com:8983", true);
  checker.checkWhitelist("abc-1.com:8983/solr", Arrays.asList(new String[]{"abc-1.com:8983/solr"}));
}

@Test
public void testWhitelistHostCheckerNonWhitelistedHostHttps() {
  WhitelistHostChecker checker = new WhitelistHostChecker("http://abc-1.com:8983, http://abc-2.com:8983, http://abc-3.com:8983", true);
  checker.checkWhitelist("https://abc-1.com:8983/solr", Arrays.asList(new String[]{"https://abc-1.com:8983/solr"}));
}

@Test
public void testWhitelistHostCheckerCoreSpecific() {
  // cores are removed completely so it doesn't really matter if they were set in config
  WhitelistHostChecker checker = new WhitelistHostChecker("http://abc-1.com:8983/solr/core1, http://abc-2.com:8983/solr2/core2", true);
  checker.checkWhitelist("http://abc-1.com:8983/solr/core2", Arrays.asList(new String[]{"http://abc-1.com:8983/solr/core2"}));
}
```

### JMeter-5.1.1 (minT = 30)

// 140	test/src	org.apache.jmeter.assertions	TestJSONPathAssertion	testGetResult_not_regexp	()V	120	135
// 140	test/src	org.apache.jmeter.assertions	TestJSONPathAssertion	testGetResult_positive_invert	()V	104	118
// 140	test/src	org.apache.jmeter.assertions	TestJSONPathAssertion	testGetResult_negative_invert	()V	152	166
// 140	test/src	org.apache.jmeter.assertions	TestJSONPathAssertion	testGetResult_inverted_null	()V	300	314
// 140	test/src	org.apache.jmeter.assertions	TestJSONPathAssertion	testGetResultFloat	()V	368	384

```java
@Test
public void testGetResult_not_regexp() {
    SampleResult samplerResult = new SampleResult();
    samplerResult.setResponseData("{\"myval\": \"some complicated value\"}".getBytes());

    JSONPathAssertion instance = new JSONPathAssertion();
    instance.setJsonPath("$.myval");
    instance.setJsonValidationBool(true);
    instance.setExpectedValue("some.+");
    AssertionResult result = instance.getResult(samplerResult);
    assertEquals(false, result.isFailure());

    instance.setIsRegex(false);
    AssertionResult result2 = instance.getResult(samplerResult);
    assertEquals(true, result2.isFailure());
}

@Test
public void testGetResult_positive_invert() {
    SampleResult samplerResult = new SampleResult();
    samplerResult.setResponseData("{\"myval\": 123}".getBytes());

    JSONPathAssertion instance = new JSONPathAssertion();
    instance.setJsonPath("$.myval");
    instance.setJsonValidationBool(true);
    instance.setExpectedValue("123");
    instance.setInvert(true);
    AssertionResult expResult = new AssertionResult("");
    AssertionResult result = instance.getResult(samplerResult);
    assertEquals(true, result.isFailure());
    assertEquals(expResult.getName(), result.getName());
}

@Test
public void testGetResult_negative_invert() {
    SampleResult samplerResult = new SampleResult();
    samplerResult.setResponseData("{\"myval\": 123}".getBytes());

    JSONPathAssertion instance = new JSONPathAssertion();
    instance.setJsonPath("$.myval");
    instance.setJsonValidationBool(true);
    instance.setExpectedValue("1234");
    instance.setInvert(true);
    AssertionResult expResult = new AssertionResult("");
    AssertionResult result = instance.getResult(samplerResult);
    assertEquals(false, result.isFailure());
    assertEquals(expResult.getName(), result.getName());
}

@Test
public void testGetResult_inverted_null() {
    SampleResult samplerResult = new SampleResult();
    samplerResult.setResponseData("{\"myval\": [{\"key\": null}]}".getBytes());

    JSONPathAssertion instance = new JSONPathAssertion();
    instance.setJsonPath("$.myval[*].key");
    instance.setJsonValidationBool(true);
    instance.setExpectNull(true);
    instance.setInvert(true);
    AssertionResult expResult = new AssertionResult("");
    AssertionResult result = instance.getResult(samplerResult);
    assertEquals(expResult.getName(), result.getName());
    assertEquals(true, result.isFailure());
}

@Test
public void testGetResultFloat() {
    SampleResult samplerResult = new SampleResult();

    samplerResult.setResponseData("{\"myval\": [{\"test\":0.0000123456789}]}".getBytes());

    JSONPathAssertion instance = new JSONPathAssertion();
    instance.setJsonPath("$.myval[*].test");
    instance.setJsonValidationBool(true);
    instance.setIsRegex(false);
    instance.setExpectedValue("0.0000123456789");

    AssertionResult expResult = new AssertionResult("");
    AssertionResult result = instance.getResult(samplerResult);
    assertEquals(expResult.getName(), result.getName());
    assertEquals(false, result.isFailure());
}
```

// 17	test/src	org.apache.jmeter.extractor	TestHtmlExtractorJSoup	testEmptyDefaultVariable	()V	96	103
// 17	test/src	org.apache.jmeter.extractor	TestHtmlExtractorJSoup	testVariableExtractionWithAttribute2	()V	140	147
// 17	test/src	org.apache.jmeter.extractor	TestRegexExtractor	testVariableExtraction0	()V	124	131
// 17	test/src	org.apache.jmeter.extractor	TestRegexExtractor	testVariableExtraction2	()V	201	208
// 17	test/src	org.apache.jmeter.extractor	TestRegexExtractor	testVariableExtraction3	()V	220	227

```java
@Test
public void testEmptyDefaultVariable() throws Exception {
    extractor.setExpression("p.missing");
    extractor.setMatchNumber(1);
    extractor.setDefaultEmptyValue(true);
    extractor.process();
    assertEquals("", vars.get("regVal"));
}

@Test
public void testVariableExtractionWithAttribute2() throws Exception {
    extractor.setExpression("a");
    extractor.setAttribute("href");
    extractor.setMatchNumber(2);
    extractor.process();
    assertEquals("http://example2.com/", vars.get("regVal"));
}

@Test
public void testVariableExtraction0() throws Exception {
    extractor.setRegex("<(value) field=\"");
    extractor.setTemplate("$1$");
    extractor.setMatchNumber(0);
    extractor.process();
    assertEquals("value", vars.get("regVal"));
}

@Test
public void testVariableExtraction2() throws Exception {
    extractor.setRegex("<value field=\"(pinposition\\d+)\">(\\d+)</value>");
    extractor.setTemplate("$1$");
    extractor.setMatchNumber(3);
    extractor.process();
    assertEquals("pinposition3", vars.get("regVal"));
}

@Test
public void testVariableExtraction3() throws Exception {
    extractor.setRegex("<value field=\"(pinposition\\d+)\">(\\d+)</value>");
    extractor.setTemplate("_$1$");
    extractor.setMatchNumber(2);
    extractor.process();
    assertEquals("_pinposition2", vars.get("regVal"));
}
```

// 198	test/src	org.apache.commons.cli.avalon	ClutilTestCase	testIncomplete2ArgsMixed	()V	813	833
// 198	test/src	org.apache.commons.cli.avalon	ClutilTestCase	testIncomplete2ArgsMixedNoEq	()V	835	855
// 198	test/src	org.apache.commons.cli.avalon	ClutilTestCase	testOptionalArgsWithArgShortBeforeOtherOpt	()V	232	254
// 198	test/src	org.apache.commons.cli.avalon	ClutilTestCase	testOptionalArgsWithArgShortEqualsBeforeOtherOpt	()V	256	278
// 198	test/src	org.apache.commons.cli.avalon	ClutilTestCase	testOptionalArgsNoArgShortBeforeOtherOpt	()V	280	302

```java
@Test
public void testIncomplete2ArgsMixed() {
    // "-Dstupid=","-c"
    final CLOptionDescriptor[] options = new CLOptionDescriptor[] { DEFINE, CLEAR1 };

    final String[] args = new String[] { "-Dstupid=", "-c" };

    final CLArgsParser parser = new CLArgsParser(args, options);

    assertNull(parser.getErrorString(), parser.getErrorString());

    final List<CLOption> clOptions = parser.getArguments();
    final int size = clOptions.size();

    assertEquals(size, 2);
    assertEquals(clOptions.get(1).getDescriptor().getId(), CLEAR1_OPT);
    final CLOption option = clOptions.get(0);
    assertEquals(option.getDescriptor().getId(), DEFINE_OPT);
    assertEquals(option.getArgument(0), "stupid");
    assertEquals(option.getArgument(1), "");
}

@Test
public void testIncomplete2ArgsMixedNoEq() {
    // "-Dstupid","-c"
    final CLOptionDescriptor[] options = new CLOptionDescriptor[] { DEFINE, CLEAR1 };

    final String[] args = new String[] { "-DStupid", "-c" };

    final CLArgsParser parser = new CLArgsParser(args, options);

    assertNull(parser.getErrorString(), parser.getErrorString());

    final List<CLOption> clOptions = parser.getArguments();
    final int size = clOptions.size();

    assertEquals(size, 2);
    assertEquals(clOptions.get(1).getDescriptor().getId(), CLEAR1_OPT);
    final CLOption option = clOptions.get(0);
    assertEquals(option.getDescriptor().getId(), DEFINE_OPT);
    assertEquals(option.getArgument(0), "Stupid");
    assertEquals(option.getArgument(1), "");
}

@Test
public void testOptionalArgsWithArgShortBeforeOtherOpt() {
    // "-T3","-a"
    final CLOptionDescriptor[] options = new CLOptionDescriptor[] { ALL, TAINT };

    final String[] args = new String[] { "-T3", "-a" };

    final CLArgsParser parser = new CLArgsParser(args, options);

    assertNull(parser.getErrorString(), parser.getErrorString());

    final List<CLOption> clOptions = parser.getArguments();
    final int size = clOptions.size();

    assertEquals(size, 2);
    final CLOption option0 = clOptions.get(0);
    assertEquals(option0.getDescriptor().getId(), TAINT_OPT);
    assertEquals(option0.getArgument(0), "3");

    final CLOption option1 = clOptions.get(1);
    assertEquals(ALL_OPT, option1.getDescriptor().getId());
    assertEquals(null, option1.getArgument(0));
}

@Test
public void testOptionalArgsWithArgShortEqualsBeforeOtherOpt() {
    // "-T3","-a"
    final CLOptionDescriptor[] options = new CLOptionDescriptor[] { ALL, TAINT };

    final String[] args = new String[] { "-T=3", "-a" };

    final CLArgsParser parser = new CLArgsParser(args, options);

    assertNull(parser.getErrorString(), parser.getErrorString());

    final List<CLOption> clOptions = parser.getArguments();
    final int size = clOptions.size();

    assertEquals(size, 2);
    final CLOption option0 = clOptions.get(0);
    assertEquals(option0.getDescriptor().getId(), TAINT_OPT);
    assertEquals(option0.getArgument(0), "3");

    final CLOption option1 = clOptions.get(1);
    assertEquals(ALL_OPT, option1.getDescriptor().getId());
    assertEquals(null, option1.getArgument(0));
}

@Test
public void testOptionalArgsNoArgShortBeforeOtherOpt() {
    // "-T","-a"
    final CLOptionDescriptor[] options = new CLOptionDescriptor[] { ALL, TAINT };

    final String[] args = new String[] { "-T", "-a" };

    final CLArgsParser parser = new CLArgsParser(args, options);

    assertNull(parser.getErrorString(), parser.getErrorString());

    final List<CLOption> clOptions = parser.getArguments();
    final int size = clOptions.size();

    assertEquals(size, 2);
    final CLOption option0 = clOptions.get(0);
    assertEquals(TAINT_OPT, option0.getDescriptor().getId());
    assertEquals(null, option0.getArgument(0));

    final CLOption option1 = clOptions.get(1);
    assertEquals(ALL_OPT, option1.getDescriptor().getId());
    assertEquals(null, option1.getArgument(0));
}
```

// 302	test/src	org.apache.jmeter.functions	TestEscapeOroRegexpChars	testEscapeWithVars	()V	88	96
// 302	test/src	org.apache.jmeter.functions	TestGroovyFunction	testSumVar	()V	80	88
// 302	test/src	org.apache.jmeter.functions	TestJavascriptFunction	testSumVar	()V	79	87
// 302	test/src	org.apache.jmeter.functions	TestJexl2Function	testSumVar	()V	72	80
// 302	test/src	org.apache.jmeter.functions	TestSetProperty	testSetPropertyNoReturn	()V	67	75

```java
@Test
public void testEscapeWithVars() throws Exception {
    params.add(new CompoundVariable("toto(.+?)titi"));
    params.add(new CompoundVariable("exportedVar"));
    function.setParameters(params);
    String ret = function.execute(result, null);
    assertEquals("toto\\(\\.\\+\\?\\)titi", ret);
    assertEquals("toto\\(\\.\\+\\?\\)titi", vars.get("exportedVar"));
}

@Test
public void testSumVar() throws Exception {
    params.add(new CompoundVariable("1+2+3"));
    params.add(new CompoundVariable("TOTAL"));
    function.setParameters(params);
    String ret = function.execute(result, null);
    assertEquals("6", ret);
    assertEquals("6", vars.get("TOTAL"));
}

@Test
public void testSumVar() throws Exception {
    params.add(new CompoundVariable("1+2+3"));
    params.add(new CompoundVariable("TOTAL"));
    function.setParameters(params);
    String ret = function.execute(result, null);
    assertEquals("6", ret);
    assertEquals("6", vars.get("TOTAL"));
}

@Test
public void testSumVar() throws Exception {
    params.add(new CompoundVariable("1+2+3"));
    params.add(new CompoundVariable("TOTAL"));
    function.setParameters(params);
    String ret = function.execute(result, null);
    assertEquals("6", ret);
    assertEquals("6", vars.get("TOTAL"));
}

@Test
public void testSetPropertyNoReturn() throws Exception {
    params.add(new CompoundVariable("prop1"));
    params.add(new CompoundVariable("value1"));
    function.setParameters(params);
    String returnValue = function.execute(result, null);
    assertEquals("value1", JMeterUtils.getProperty("prop1"));
    assertEquals("", returnValue);
}
```

### commons-collection-4.3 (minT = 30)

// also included in minT = 100

// 323	src/test/java	org.apache.commons.collections4.map	Flat3MapTest	testMapIteratorSetValue1	()V	279	296
// 323	src/test/java	org.apache.commons.collections4.map	Flat3MapTest	testMapIteratorSetValue2	()V	298	316
// 323	src/test/java	org.apache.commons.collections4.map	Flat3MapTest	testEntryIteratorSetValue1	()V	218	235
// 323	src/test/java	org.apache.commons.collections4.map	Flat3MapTest	testEntryIteratorSetValue2	()V	237	255
// 323	src/test/java	org.apache.commons.collections4.map	Flat3MapTest	testEntryIteratorSetValue3	()V	257	276

```java
@SuppressWarnings("unchecked")
public void testMapIteratorSetValue1() throws Exception {
    final Flat3Map<K, V> map = makeObject();
    map.put((K) ONE, (V) TEN);
    map.put((K) TWO, (V) TWENTY);
    map.put((K) THREE, (V) THIRTY);

    final MapIterator<K, V> it = map.mapIterator();
    it.next();
    it.setValue((V) "NewValue");
    assertEquals(3, map.size());
    assertEquals(true, map.containsKey(ONE));
    assertEquals(true, map.containsKey(TWO));
    assertEquals(true, map.containsKey(THREE));
    assertEquals("NewValue", map.get(ONE));
    assertEquals(TWENTY, map.get(TWO));
    assertEquals(THIRTY, map.get(THREE));
}

@SuppressWarnings("unchecked")
public void testMapIteratorSetValue2() throws Exception {
    final Flat3Map<K, V> map = makeObject();
    map.put((K) ONE, (V) TEN);
    map.put((K) TWO, (V) TWENTY);
    map.put((K) THREE, (V) THIRTY);

    final MapIterator<K, V> it = map.mapIterator();
    it.next();
    it.next();
    it.setValue((V) "NewValue");
    assertEquals(3, map.size());
    assertEquals(true, map.containsKey(ONE));
    assertEquals(true, map.containsKey(TWO));
    assertEquals(true, map.containsKey(THREE));
    assertEquals(TEN, map.get(ONE));
    assertEquals("NewValue", map.get(TWO));
    assertEquals(THIRTY, map.get(THREE));
}

@SuppressWarnings("unchecked")
public void testEntryIteratorSetValue1() throws Exception {
    final Flat3Map<K, V> map = makeObject();
    map.put((K) ONE, (V) TEN);
    map.put((K) TWO, (V) TWENTY);
    map.put((K) THREE, (V) THIRTY);

    final Iterator<Map.Entry<K, V>> it = map.entrySet().iterator();
    final Map.Entry<K, V> entry = it.next();
    entry.setValue((V) "NewValue");
    assertEquals(3, map.size());
    assertEquals(true, map.containsKey(ONE));
    assertEquals(true, map.containsKey(TWO));
    assertEquals(true, map.containsKey(THREE));
    assertEquals("NewValue", map.get(ONE));
    assertEquals(TWENTY, map.get(TWO));
    assertEquals(THIRTY, map.get(THREE));
}

@SuppressWarnings("unchecked")
public void testEntryIteratorSetValue2() throws Exception {
    final Flat3Map<K, V> map = makeObject();
    map.put((K) ONE, (V) TEN);
    map.put((K) TWO, (V) TWENTY);
    map.put((K) THREE, (V) THIRTY);

    final Iterator<Map.Entry<K, V>> it = map.entrySet().iterator();
    it.next();
    final Map.Entry<K, V> entry = it.next();
    entry.setValue((V) "NewValue");
    assertEquals(3, map.size());
    assertEquals(true, map.containsKey(ONE));
    assertEquals(true, map.containsKey(TWO));
    assertEquals(true, map.containsKey(THREE));
    assertEquals(TEN, map.get(ONE));
    assertEquals("NewValue", map.get(TWO));
    assertEquals(THIRTY, map.get(THREE));
}

@SuppressWarnings("unchecked")
public void testEntryIteratorSetValue3() throws Exception {
    final Flat3Map<K, V> map = makeObject();
    map.put((K) ONE, (V) TEN);
    map.put((K) TWO, (V) TWENTY);
    map.put((K) THREE, (V) THIRTY);

    final Iterator<Map.Entry<K, V>> it = map.entrySet().iterator();
    it.next();
    it.next();
    final Map.Entry<K, V> entry = it.next();
    entry.setValue((V) "NewValue");
    assertEquals(3, map.size());
    assertEquals(true, map.containsKey(ONE));
    assertEquals(true, map.containsKey(TWO));
    assertEquals(true, map.containsKey(THREE));
    assertEquals(TEN, map.get(ONE));
    assertEquals(TWENTY, map.get(TWO));
    assertEquals("NewValue", map.get(THREE));
}
```

### commons-io-2.5

// 136	src/test/java	org.apache.commons.io	FilenameUtilsTestCase	testGetPrefix_with_nullbyte	()V	647	654
// 136	src/test/java	org.apache.commons.io	FilenameUtilsTestCase	testGetPathNoEndSeparator_with_null_byte	()V	736	743
// 136	src/test/java	org.apache.commons.io	FilenameUtilsTestCase	testInjectionFailure	()V	844	851
// 136	src/test/java	org.apache.commons.io	FilenameUtilsTestCase	testGetBaseName_with_nullByte	()V	864	871
// 136	src/test/java	org.apache.commons.io	FilenameUtilsTestCase	testIsExtension_injection	()V	1022	1029

```java
@Test
public void testGetPrefix_with_nullbyte() {
    try {
        assertEquals("~user\\", FilenameUtils.getPrefix("~u\u0000ser\\a\\b\\c.txt"));
    } catch (IllegalArgumentException ignore) {

    }
}

@Test
public void testGetPathNoEndSeparator_with_null_byte() {
    try {
        assertEquals("a/b", FilenameUtils.getPathNoEndSeparator("~user/a\u0000/b/c.txt"));
    } catch (IllegalArgumentException ignore) {

    }
}

@Test
public void testInjectionFailure() {
    try {
        assertEquals("c", FilenameUtils.getName("a\\b\\\u0000c"));
    } catch (IllegalArgumentException ignore) {

    }
}

@Test
public void testGetBaseName_with_nullByte() {
    try {
        assertEquals("file.txt", FilenameUtils.getBaseName("fil\u0000e.txt.bak"));
    } catch (IllegalArgumentException ignore) {

    }
}

@Test
public void testIsExtension_injection() {
    try {
        FilenameUtils.isExtension("a.b\\fi\u0000le.txt", "TXT");
        fail("Should throw IAE");
    } catch (IllegalArgumentException ignore) {
    }
}
```

// 142	src/test/java	org.apache.commons.io	CopyUtilsTest	copy_byteArrayToWriter	()V	82	93
// 142	src/test/java	org.apache.commons.io	IOUtilsWriteTestCase	testWrite_stringToOutputStream	()V	317	330
// 142	src/test/java	org.apache.commons.io	IOUtilsWriteTestCase	testWrite_charArrayToOutputStream	()V	451	464
// 142	src/test/java	org.apache.commons.io	IOUtilsWriteTestCase	testWrite_byteArrayToWriter	()V	89	102
// 142	src/test/java	org.apache.commons.io	IOUtilsWriteTestCase	testWrite_charSequenceToOutputStream	()V	183	196

```java
@Test
public void copy_byteArrayToWriter() throws Exception {
    final ByteArrayOutputStream baout = new ByteArrayOutputStream();
    final OutputStream out = new YellOnFlushAndCloseOutputStream(baout, false, true);
    final Writer writer = new java.io.OutputStreamWriter(out, "US-ASCII");

    CopyUtils.copy(inData, writer);
    writer.flush();

    assertEquals("Sizes differ", inData.length, baout.size());
    assertTrue("Content differs", Arrays.equals(inData, baout.toByteArray()));
}

@Test
public void testWrite_stringToOutputStream() throws Exception {
    final String str = new String(inData, "US-ASCII");

    final ByteArrayOutputStream baout = new ByteArrayOutputStream();
    final YellOnFlushAndCloseOutputStream out = new YellOnFlushAndCloseOutputStream(baout, true, true);

    IOUtils.write(str, out);
    out.off();
    out.flush();

    assertEquals("Sizes differ", inData.length, baout.size());
    assertTrue("Content differs", Arrays.equals(inData, baout.toByteArray()));
}

@Test
public void testWrite_charArrayToOutputStream() throws Exception {
    final String str = new String(inData, "US-ASCII");

    final ByteArrayOutputStream baout = new ByteArrayOutputStream();
    final YellOnFlushAndCloseOutputStream out = new YellOnFlushAndCloseOutputStream(baout, true, true);

    IOUtils.write(str.toCharArray(), out);
    out.off();
    out.flush();

    assertEquals("Sizes differ", inData.length, baout.size());
    assertTrue("Content differs", Arrays.equals(inData, baout.toByteArray()));
}

@Test
public void testWrite_byteArrayToWriter() throws Exception {
    final ByteArrayOutputStream baout = new ByteArrayOutputStream();
    @SuppressWarnings("resource") // deliberately not closed
    final YellOnFlushAndCloseOutputStream out = new YellOnFlushAndCloseOutputStream(baout, true, true);
    final Writer writer = new OutputStreamWriter(baout, "US-ASCII");

    IOUtils.write(inData, writer);
    out.off();
    writer.flush();

    assertEquals("Sizes differ", inData.length, baout.size());
    assertTrue("Content differs", Arrays.equals(inData, baout.toByteArray()));
}

@Test
public void testWrite_charSequenceToOutputStream() throws Exception {
    final CharSequence csq = new StringBuilder(new String(inData, "US-ASCII"));

    final ByteArrayOutputStream baout = new ByteArrayOutputStream();
    final YellOnFlushAndCloseOutputStream out = new YellOnFlushAndCloseOutputStream(baout, true, true);

    IOUtils.write(csq, out);
    out.off();
    out.flush();

    assertEquals("Sizes differ", inData.length, baout.size());
    assertTrue("Content differs", Arrays.equals(inData, baout.toByteArray()));
}
```

// 232	src/test/java	org.apache.commons.io	FileSystemUtilsTestCase	testGetFreeSpaceUnix_String_NormalResponseLinux	()V	318	326
// 232	src/test/java	org.apache.commons.io	FileSystemUtilsTestCase	testGetFreeSpaceUnix_String_NormalResponseFreeBSD	()V	328	336
// 232	src/test/java	org.apache.commons.io	FileSystemUtilsTestCase	testGetFreeSpaceUnix_String_NormalResponseKbLinux	()V	339	348
// 232	src/test/java	org.apache.commons.io	FileSystemUtilsTestCase	testGetFreeSpaceUnix_String_NormalResponseKbFreeBSD	()V	350	359
// 232	src/test/java	org.apache.commons.io	FileSystemUtilsTestCase	testGetFreeSpaceUnix_String_NormalResponseKbSolaris	()V	361	370

```java
@Test
public void testGetFreeSpaceUnix_String_NormalResponseLinux() throws Exception {
    // from Sourceforge 'GNU bash, version 2.05b.0(1)-release (i386-redhat-linux-gnu)'
    final String lines =
            "Filesystem           1K-blocks      Used Available Use% Mounted on\n" +
                    "/dev/xxx                497944    308528    189416  62% /";
    final FileSystemUtils fsu = new MockFileSystemUtils(0, lines);
    assertEquals(189416L, fsu.freeSpaceUnix("/", false, false, -1));
}

@Test
public void testGetFreeSpaceUnix_String_NormalResponseFreeBSD() throws Exception {
    // from Apache 'FreeBSD 6.1-RELEASE (SMP-turbo)'
    final String lines =
            "Filesystem  1K-blocks      Used    Avail Capacity  Mounted on\n" +
                    "/dev/xxxxxx    128990    102902    15770    87%    /";
    final FileSystemUtils fsu = new MockFileSystemUtils(0, lines);
    assertEquals(15770L, fsu.freeSpaceUnix("/", false, false, -1));
}

//-----------------------------------------------------------------------
@Test
public void testGetFreeSpaceUnix_String_NormalResponseKbLinux() throws Exception {
    // from Sourceforge 'GNU bash, version 2.05b.0(1)-release (i386-redhat-linux-gnu)'
    // df, df -k and df -kP are all identical
    final String lines =
            "Filesystem           1K-blocks      Used Available Use% Mounted on\n" +
                    "/dev/xxx                497944    308528    189416  62% /";
    final FileSystemUtils fsu = new MockFileSystemUtils(0, lines);
    assertEquals(189416L, fsu.freeSpaceUnix("/", true, false, -1));
}

@Test
public void testGetFreeSpaceUnix_String_NormalResponseKbFreeBSD() throws Exception {
    // from Apache 'FreeBSD 6.1-RELEASE (SMP-turbo)'
    // df and df -k are identical, but df -kP uses 512 blocks (not relevant as not used)
    final String lines =
            "Filesystem  1K-blocks      Used    Avail Capacity  Mounted on\n" +
                    "/dev/xxxxxx    128990    102902    15770    87%    /";
    final FileSystemUtils fsu = new MockFileSystemUtils(0, lines);
    assertEquals(15770L, fsu.freeSpaceUnix("/", true, false, -1));
}

@Test
public void testGetFreeSpaceUnix_String_NormalResponseKbSolaris() throws Exception {
    // from IO-91 - ' SunOS et 5.10 Generic_118822-25 sun4u sparc SUNW,Ultra-4'
    // non-kb response does not contain free space - see IO-91
    final String lines =
            "Filesystem            kbytes    used   avail capacity  Mounted on\n" +
                    "/dev/dsk/x0x0x0x0    1350955  815754  481163    63%";
    final FileSystemUtils fsu = new MockFileSystemUtils(0, lines);
    assertEquals(481163L, fsu.freeSpaceUnix("/dev/dsk/x0x0x0x0", true, false, -1));
}
```

// 430	src/test/java	org.apache.commons.io	FileUtilsTestCase	testWriteStringToFile1	()V	1979	1985
// 430	src/test/java	org.apache.commons.io	FileUtilsTestCase	testWriteCharSequence1	()V	2003	2009
// 430	src/test/java	org.apache.commons.io	FileUtilsTestCase	testWriteStringToFile2	()V	1987	1993
// 430	src/test/java	org.apache.commons.io	FileUtilsTestCase	testWriteStringToFile3	()V	1995	2001
// 430	src/test/java	org.apache.commons.io	FileUtilsTestCase	testWriteCharSequence2	()V	2011	2017

```java
@Test
public void testWriteStringToFile1() throws Exception {
    final File file = new File(getTestDirectory(), "write.txt");
    FileUtils.writeStringToFile(file, "Hello /u1234", "UTF8");
    final byte[] text = "Hello /u1234".getBytes("UTF8");
    TestUtils.assertEqualContent(text, file);
}

@Test
public void testWriteCharSequence1() throws Exception {
    final File file = new File(getTestDirectory(), "write.txt");
    FileUtils.write(file, "Hello /u1234", "UTF8");
    final byte[] text = "Hello /u1234".getBytes("UTF8");
    TestUtils.assertEqualContent(text, file);
}

@Test
public void testWriteStringToFile2() throws Exception {
    final File file = new File(getTestDirectory(), "write.txt");
    FileUtils.writeStringToFile(file, "Hello /u1234", (String) null);
    final byte[] text = "Hello /u1234".getBytes();
    TestUtils.assertEqualContent(text, file);
}

@Test
public void testWriteStringToFile3() throws Exception {
    final File file = new File(getTestDirectory(), "write.txt");
    FileUtils.writeStringToFile(file, "Hello /u1234", (Charset) null);
    final byte[] text = "Hello /u1234".getBytes();
    TestUtils.assertEqualContent(text, file);
}

@Test
public void testWriteCharSequence2() throws Exception {
    final File file = new File(getTestDirectory(), "write.txt");
    FileUtils.write(file, "Hello /u1234", (String) null);
    final byte[] text = "Hello /u1234".getBytes();
    TestUtils.assertEqualContent(text, file);
}
```