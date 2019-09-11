[8, 'src/test/java', 'org.apache.commons.math3.linear', 'BlockFieldMatrixTest', 'testGetVectors', 482, 499]

[8, 'src/test/java', 'org.apache.commons.math3.linear', 'FieldMatrixImplTest', 'testGetVectors', 345, 362]

    @Test
    public void testGetVectors() {
        FieldMatrix<Fraction> m = new BlockFieldMatrix<Fraction>(testData);
        TestUtils.assertEquals(m.getRow(0), testDataRow1);
        TestUtils.assertEquals(m.getColumn(2), testDataCol3);
        try {
            m.getRow(10);
            Assert.fail("expecting OutOfRangeException");
        } catch (OutOfRangeException ex) {
            // ignored
        }
        try {
            m.getColumn(-1);
            Assert.fail("expecting OutOfRangeException");
        } catch (OutOfRangeException ex) {
            // ignored
        }
    }

    @Test
    public void testGetVectors() {
        FieldMatrix<Fraction> m = new Array2DRowFieldMatrix<Fraction>(testData);
        TestUtils.assertEquals(m.getRow(0), testDataRow1);
        TestUtils.assertEquals(m.getColumn(2), testDataCol3);
        try {
            m.getRow(10);
            Assert.fail("expecting OutOfRangeException");
        } catch (OutOfRangeException ex) {
            // ignored
        }
        try {
            m.getColumn(-1);
            Assert.fail("expecting OutOfRangeException");
        } catch (OutOfRangeException ex) {
            // ignored
        }
    }

[19, 'src/test/java', 'org.apache.commons.math3.linear', 'BlockFieldMatrixTest', 'testMultiply', 226, 243]

[19, 'src/test/java', 'org.apache.commons.math3.linear', 'FieldMatrixImplTest', 'testMultiply', 167, 184]

    @Test
    public void testMultiply() {
        BlockFieldMatrix<Fraction> m = new BlockFieldMatrix<Fraction>(testData);
        BlockFieldMatrix<Fraction> mInv = new BlockFieldMatrix<Fraction>(testDataInv);
        BlockFieldMatrix<Fraction> identity = new BlockFieldMatrix<Fraction>(id);
        BlockFieldMatrix<Fraction> m2 = new BlockFieldMatrix<Fraction>(testData2);
        TestUtils.assertEquals(m.multiply(mInv), identity);
        TestUtils.assertEquals(mInv.multiply(m), identity);
        TestUtils.assertEquals(m.multiply(identity), m);
        TestUtils.assertEquals(identity.multiply(mInv), mInv);
        TestUtils.assertEquals(m2.multiply(identity), m2);
        try {
            m.multiply(new BlockFieldMatrix<Fraction>(bigSingular));
            Assert.fail("Expecting illegalArgumentException");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }
    }

    @Test
     public void testMultiply() {
        Array2DRowFieldMatrix<Fraction> m = new Array2DRowFieldMatrix<Fraction>(testData);
        Array2DRowFieldMatrix<Fraction> mInv = new Array2DRowFieldMatrix<Fraction>(testDataInv);
        Array2DRowFieldMatrix<Fraction> identity = new Array2DRowFieldMatrix<Fraction>(id);
        Array2DRowFieldMatrix<Fraction> m2 = new Array2DRowFieldMatrix<Fraction>(testData2);
        TestUtils.assertEquals(m.multiply(mInv), identity);
        TestUtils.assertEquals(mInv.multiply(m), identity);
        TestUtils.assertEquals(m.multiply(identity), m);
        TestUtils.assertEquals(identity.multiply(mInv), mInv);
        TestUtils.assertEquals(m2.multiply(identity), m2);
        try {
            m.multiply(new Array2DRowFieldMatrix<Fraction>(bigSingular));
            Assert.fail("Expecting illegalArgumentException");
        } catch (MathIllegalArgumentException ex) {
            // ignored
        }
    }

[20, 'src/test/java', 'org.apache.commons.math3.linear', 'BlockRealMatrixTest', 'testGetSetRowLarge', 952, 970]

[20, 'src/test/java', 'org.apache.commons.math3.linear', 'BlockRealMatrixTest', 'testGetSetColumnLarge', 1014, 1032]

    @Test
    public void testGetSetRowLarge() {
        int n = 3 * BlockRealMatrix.BLOCK_SIZE;
        RealMatrix m = new BlockRealMatrix(n, n);
        double[] sub = new double[n];
        Arrays.fill(sub, 1.0);

        m.setRow(2, sub);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                if (i != 2) {
                    Assert.assertEquals(0.0, m.getEntry(i, j), 0.0);
                } else {
                    Assert.assertEquals(1.0, m.getEntry(i, j), 0.0);
                }
            }
        }
        checkArrays(sub, m.getRow(2));
    }

    @Test
    public void testGetSetColumnLarge() {
        int n = 3 * BlockRealMatrix.BLOCK_SIZE;
        RealMatrix m = new BlockRealMatrix(n, n);
        double[] sub = new double[n];
        Arrays.fill(sub, 1.0);

        m.setColumn(2, sub);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                if (j != 2) {
                    Assert.assertEquals(0.0, m.getEntry(i, j), 0.0);
                } else {
                    Assert.assertEquals(1.0, m.getEntry(i, j), 0.0);
                }
            }
        }
        checkArrays(sub, m.getColumn(2));
    }

[54, 'src/test/java', 'org.apache.commons.math3.stat.inference', 'GTestTest', 'testRootLogLikelihood', 270, 290]

[54, 'src/test/java', 'org.apache.commons.math3.stat.inference', 'TestUtilsTest', 'testRootLogLikelihood', 510, 530]

    @Test
    public void testRootLogLikelihood() {
        // positive where k11 is bigger than expected.
        Assert.assertTrue(testStatistic.rootLogLikelihoodRatio(904, 21060, 1144, 283012) > 0.0);

        // negative because k11 is lower than expected
        Assert.assertTrue(testStatistic.rootLogLikelihoodRatio(36, 21928, 60280, 623876) < 0.0);

        Assert.assertEquals(FastMath.sqrt(2.772589), testStatistic.rootLogLikelihoodRatio(1, 0, 0, 1), 0.000001);
        Assert.assertEquals(-FastMath.sqrt(2.772589), testStatistic.rootLogLikelihoodRatio(0, 1, 1, 0), 0.000001);
        Assert.assertEquals(FastMath.sqrt(27.72589), testStatistic.rootLogLikelihoodRatio(10, 0, 0, 10), 0.00001);

        Assert.assertEquals(FastMath.sqrt(39.33052), testStatistic.rootLogLikelihoodRatio(5, 1995, 0, 100000), 0.00001);
        Assert.assertEquals(-FastMath.sqrt(39.33052), testStatistic.rootLogLikelihoodRatio(0, 100000, 5, 1995), 0.00001);

        Assert.assertEquals(FastMath.sqrt(4730.737), testStatistic.rootLogLikelihoodRatio(1000, 1995, 1000, 100000), 0.001);
        Assert.assertEquals(-FastMath.sqrt(4730.737), testStatistic.rootLogLikelihoodRatio(1000, 100000, 1000, 1995), 0.001);

        Assert.assertEquals(FastMath.sqrt(5734.343), testStatistic.rootLogLikelihoodRatio(1000, 1000, 1000, 100000), 0.001);
        Assert.assertEquals(FastMath.sqrt(5714.932), testStatistic.rootLogLikelihoodRatio(1000, 1000, 1000, 99000), 0.001);
    }

    @Test
    public void testRootLogLikelihood() {
        // positive where k11 is bigger than expected.
        Assert.assertTrue(TestUtils.rootLogLikelihoodRatio(904, 21060, 1144, 283012) > 0.0);

        // negative because k11 is lower than expected
        Assert.assertTrue(TestUtils.rootLogLikelihoodRatio(36, 21928, 60280, 623876) < 0.0);

        Assert.assertEquals(FastMath.sqrt(2.772589), TestUtils.rootLogLikelihoodRatio(1, 0, 0, 1), 0.000001);
        Assert.assertEquals(-FastMath.sqrt(2.772589), TestUtils.rootLogLikelihoodRatio(0, 1, 1, 0), 0.000001);
        Assert.assertEquals(FastMath.sqrt(27.72589), TestUtils.rootLogLikelihoodRatio(10, 0, 0, 10), 0.00001);

        Assert.assertEquals(FastMath.sqrt(39.33052), TestUtils.rootLogLikelihoodRatio(5, 1995, 0, 100000), 0.00001);
        Assert.assertEquals(-FastMath.sqrt(39.33052), TestUtils.rootLogLikelihoodRatio(0, 100000, 5, 1995), 0.00001);

        Assert.assertEquals(FastMath.sqrt(4730.737), TestUtils.rootLogLikelihoodRatio(1000, 1995, 1000, 100000), 0.001);
        Assert.assertEquals(-FastMath.sqrt(4730.737), TestUtils.rootLogLikelihoodRatio(1000, 100000, 1000, 1995), 0.001);

        Assert.assertEquals(FastMath.sqrt(5734.343), TestUtils.rootLogLikelihoodRatio(1000, 1000, 1000, 100000), 0.001);
        Assert.assertEquals(FastMath.sqrt(5714.932), TestUtils.rootLogLikelihoodRatio(1000, 1000, 1000, 99000), 0.001);
    }

[73, 'src/test/java', 'org.apache.commons.math3.random', 'RandomDataGeneratorTest', 'testNextIntNegativeRange', 116, 123]

[73, 'src/test/java', 'org.apache.commons.math3.random', 'RandomDataGeneratorTest', 'testNextLongNegativeRange', 191, 198]

    @Test
    public void testNextIntNegativeRange() {
        for (int i = 0; i < 5; i++) {
            checkNextIntUniform(-7, -4);
            checkNextIntUniform(-15, -2);
            checkNextIntUniform(Integer.MIN_VALUE + 1, Integer.MIN_VALUE + 12);
        }
    }

    @Test
    public void testNextLongNegativeRange() {
        for (int i = 0; i < 5; i++) {
            checkNextLongUniform(-7, -4);
            checkNextLongUniform(-15, -2);
            checkNextLongUniform(Long.MIN_VALUE + 1, Long.MIN_VALUE + 12);
        }
    }

[84, 'src/test/java', 'org.apache.commons.math3.linear', 'MatrixUtilsTest', 'testBigFractionConverter', 252, 263]

[84, 'src/test/java', 'org.apache.commons.math3.linear', 'MatrixUtilsTest', 'testFractionConverter', 265, 276]

    @Test
    public void testBigFractionConverter() {
        BigFraction[][] bfData = {
                { new BigFraction(1), new BigFraction(2), new BigFraction(3) },
                { new BigFraction(2), new BigFraction(5), new BigFraction(3) },
                { new BigFraction(1), new BigFraction(0), new BigFraction(8) }
        };
        FieldMatrix<BigFraction> m = new Array2DRowFieldMatrix<BigFraction>(bfData, false);
        RealMatrix converted = MatrixUtils.bigFractionMatrixToRealMatrix(m);
        RealMatrix reference = new Array2DRowRealMatrix(testData, false);
        Assert.assertEquals(0.0, converted.subtract(reference).getNorm(), 0.0);
    }

    @Test
    public void testFractionConverter() {
        Fraction[][] fData = {
                { new Fraction(1), new Fraction(2), new Fraction(3) },
                { new Fraction(2), new Fraction(5), new Fraction(3) },
                { new Fraction(1), new Fraction(0), new Fraction(8) }
        };
        FieldMatrix<Fraction> m = new Array2DRowFieldMatrix<Fraction>(fData, false);
        RealMatrix converted = MatrixUtils.fractionMatrixToRealMatrix(m);
        RealMatrix reference = new Array2DRowRealMatrix(testData, false);
        Assert.assertEquals(0.0, converted.subtract(reference).getNorm(), 0.0);
    }

[89, 'src/test/java', 'org.apache.commons.math3.ode.nonstiff', 'DormandPrince853IntegratorTest', 'testUnstableDerivative', 323, 334]

[89, 'src/test/java', 'org.apache.commons.math3.ode.nonstiff', 'GraggBulirschStoerIntegratorTest', 'testUnstableDerivative', 301, 312]

  @Test
  public void testUnstableDerivative()
      throws DimensionMismatchException, NumberIsTooSmallException,
             MaxCountExceededException, NoBracketingException {
    final StepProblem stepProblem = new StepProblem(0.0, 1.0, 2.0);
    FirstOrderIntegrator integ =
      new DormandPrince853Integrator(0.1, 10, 1.0e-12, 0.0);
    integ.addEventHandler(stepProblem, 1.0, 1.0e-12, 1000);
    double[] y = { Double.NaN };
    integ.integrate(stepProblem, 0.0, new double[] { 0.0 }, 10.0, y);
    Assert.assertEquals(8.0, y[0], 1.0e-12);
  }

  @Test
  public void testUnstableDerivative()
      throws DimensionMismatchException, NumberIsTooSmallException,
             MaxCountExceededException, NoBracketingException {
    final StepProblem stepProblem = new StepProblem(0.0, 1.0, 2.0);
    FirstOrderIntegrator integ =
      new GraggBulirschStoerIntegrator(0.1, 10, 1.0e-12, 0.0);
    integ.addEventHandler(stepProblem, 1.0, 1.0e-12, 1000);
    double[] y = { Double.NaN };
    integ.integrate(stepProblem, 0.0, new double[] { 0.0 }, 10.0, y);
    Assert.assertEquals(8.0, y[0], 1.0e-12);
  }

[126, 'src/test/java', 'org.apache.commons.math3.linear', 'Array2DRowRealMatrixTest', 'testSetColumnVector', 765, 784]

[126, 'src/test/java', 'org.apache.commons.math3.linear', 'BlockRealMatrixTest', 'testSetColumnVector', 865, 884]

    @Test
    public void testSetColumnVector() {
        RealMatrix m = new Array2DRowRealMatrix(subTestData);
        RealVector mColumn3 = columnToVector(subColumn3);
        Assert.assertNotSame(mColumn3, m.getColumnVector(1));
        m.setColumnVector(1, mColumn3);
        Assert.assertEquals(mColumn3, m.getColumnVector(1));
        try {
            m.setColumnVector(-1, mColumn3);
            Assert.fail("Expecting OutOfRangeException");
        } catch (OutOfRangeException ex) {
            // expected
        }
        try {
            m.setColumnVector(0, new ArrayRealVector(5));
            Assert.fail("Expecting MatrixDimensionMismatchException");
        } catch (MatrixDimensionMismatchException ex) {
            // expected
        }
    }

    @Test
    public void testSetColumnVector() {
        RealMatrix m = new BlockRealMatrix(subTestData);
        RealVector mColumn3 = columnToVector(subColumn3);
        Assert.assertNotSame(mColumn3, m.getColumnVector(1));
        m.setColumnVector(1, mColumn3);
        Assert.assertEquals(mColumn3, m.getColumnVector(1));
        try {
            m.setColumnVector(-1, mColumn3);
            Assert.fail("Expecting OutOfRangeException");
        } catch (OutOfRangeException ex) {
            // expected
        }
        try {
            m.setColumnVector(0, new ArrayRealVector(5));
            Assert.fail("Expecting MatrixDimensionMismatchException");
        } catch (MatrixDimensionMismatchException ex) {
            // expected
        }
    }

[152, 'src/test/java', 'org.apache.commons.math3.linear', 'BlockFieldMatrixTest', 'testGetColumn', 1061, 1080]

[152, 'src/test/java', 'org.apache.commons.math3.linear', 'FieldMatrixImplTest', 'testGetColumn', 800, 819]

    @Test
    public void testGetColumn() {
        FieldMatrix<Fraction> m = new BlockFieldMatrix<Fraction>(subTestData);
        Fraction[] mColumn1 = columnToArray(subColumn1);
        Fraction[] mColumn3 = columnToArray(subColumn3);
        checkArrays(mColumn1, m.getColumn(1));
        checkArrays(mColumn3, m.getColumn(3));
        try {
            m.getColumn(-1);
            Assert.fail("Expecting OutOfRangeException");
        } catch (OutOfRangeException ex) {
            // expected
        }
        try {
            m.getColumn(4);
            Assert.fail("Expecting OutOfRangeException");
        } catch (OutOfRangeException ex) {
            // expected
        }
    }

    @Test
    public void testGetColumn() {
        FieldMatrix<Fraction> m = new Array2DRowFieldMatrix<Fraction>(subTestData);
        Fraction[] mColumn1 = columnToArray(subColumn1);
        Fraction[] mColumn3 = columnToArray(subColumn3);
        checkArrays(mColumn1, m.getColumn(1));
        checkArrays(mColumn3, m.getColumn(3));
        try {
            m.getColumn(-1);
            Assert.fail("Expecting OutOfRangeException");
        } catch (OutOfRangeException ex) {
            // expected
        }
        try {
            m.getColumn(4);
            Assert.fail("Expecting OutOfRangeException");
        } catch (OutOfRangeException ex) {
            // expected
        }
    }

[154, 'src/test/java', 'org.apache.commons.math3.genetics', 'NPointCrossoverTest', 'testCrossoverInvalidFixedLengthChromosomeFirst', 53, 66]

[154, 'src/test/java', 'org.apache.commons.math3.genetics', 'NPointCrossoverTest', 'testCrossoverInvalidFixedLengthChromosomeSecond', 68, 81]

    @Test(expected = MathIllegalArgumentException.class)
    public void testCrossoverInvalidFixedLengthChromosomeFirst() {
        final Integer[] p1 = new Integer[] {1,0,1,0,0,1,0,1,1};
        final BinaryChromosome p1c = new DummyBinaryChromosome(p1);
        final Chromosome p2c = new Chromosome() {
            public double fitness() {
                // Not important
                return 0;
            }
        };

        final CrossoverPolicy cp = new NPointCrossover<Integer>(1);
        cp.crossover(p1c,p2c);
    }

    @Test(expected = MathIllegalArgumentException.class)
    public void testCrossoverInvalidFixedLengthChromosomeSecond() {
        final Integer[] p1 = new Integer[] {1,0,1,0,0,1,0,1,1};
        final BinaryChromosome p2c = new DummyBinaryChromosome(p1);
        final Chromosome p1c = new Chromosome() {
            public double fitness() {
                // Not important
                return 0;
            }
        };

        final CrossoverPolicy cp = new NPointCrossover<Integer>(1);
        cp.crossover(p1c,p2c);
    }

[162, 'src/test/java', 'org.apache.commons.math3.geometry.euclidean.threed', 'RotationTest', 'testApplyTo', 552, 568]

[162, 'src/test/java', 'org.apache.commons.math3.geometry.euclidean.threed', 'RotationTest', 'testApplyInverseToRotation', 608, 624]

  @Test
  public void testApplyTo() throws MathIllegalArgumentException {

    Rotation r1 = new Rotation(new Vector3D(2, -3, 5), 1.7, RotationConvention.VECTOR_OPERATOR);
    Rotation r2 = new Rotation(new Vector3D(-1, 3, 2), 0.3, RotationConvention.VECTOR_OPERATOR);
    Rotation r3 = r2.applyTo(r1);

    for (double x = -0.9; x < 0.9; x += 0.2) {
      for (double y = -0.9; y < 0.9; y += 0.2) {
        for (double z = -0.9; z < 0.9; z += 0.2) {
          Vector3D u = new Vector3D(x, y, z);
          checkVector(r2.applyTo(r1.applyTo(u)), r3.applyTo(u));
        }
      }
    }

  }

  @Test
  public void testApplyInverseToRotation() throws MathIllegalArgumentException {

    Rotation r1 = new Rotation(new Vector3D(2, -3, 5), 1.7, RotationConvention.VECTOR_OPERATOR);
    Rotation r2 = new Rotation(new Vector3D(-1, 3, 2), 0.3, RotationConvention.VECTOR_OPERATOR);
    Rotation r3 = r2.applyInverseTo(r1);

    for (double x = -0.9; x < 0.9; x += 0.2) {
      for (double y = -0.9; y < 0.9; y += 0.2) {
        for (double z = -0.9; z < 0.9; z += 0.2) {
          Vector3D u = new Vector3D(x, y, z);
          checkVector(r2.applyInverseTo(r1.applyTo(u)), r3.applyTo(u));
        }
      }
    }

  }

[175, 'src/test/java', 'org.apache.commons.math3.analysis.differentiation', 'SparseGradientTest', 'testExpm1Definition', 531, 540]

[175, 'src/test/java', 'org.apache.commons.math3.analysis.differentiation', 'SparseGradientTest', 'testLog1pDefinition', 550, 559]

    @Test
    public void testExpm1Definition() {
        for (double x = 0.1; x < 1.2; x += 0.001) {
            SparseGradient sgX = SparseGradient.createVariable(0, x);
            SparseGradient expm11 = sgX.expm1();
            SparseGradient expm12 = sgX.exp().subtract(sgX.getField().getOne());
            SparseGradient zero = expm11.subtract(expm12);
            checkF0F1(zero, 0.0, 0.0);
        }
    }

    @Test
    public void testLog1pDefinition() {
        for (double x = 0.1; x < 1.2; x += 0.001) {
            SparseGradient sgX = SparseGradient.createVariable(0, x);
            SparseGradient log1p1 = sgX.log1p();
            SparseGradient log1p2 = sgX.add(sgX.getField().getOne()).log();
            SparseGradient zero = log1p1.subtract(log1p2);
            checkF0F1(zero, 0.0, 0.0);
        }
    }

[185, 'src/test/java', 'org.apache.commons.math3.geometry.euclidean.threed', 'FieldRotationDSTest', 'testIssue801', 1202, 1219]

[185, 'src/test/java', 'org.apache.commons.math3.geometry.euclidean.threed', 'FieldRotationDfpTest', 'testIssue801', 967, 984]

    @Test
    public void testIssue801() throws MathArithmeticException {
        FieldVector3D<DerivativeStructure> u1 = createVector(0.9999988431610581, -0.0015210774290851095, 0.0);
        FieldVector3D<DerivativeStructure> u2 = createVector(0.0, 0.0, 1.0);

        FieldVector3D<DerivativeStructure> v1 = createVector(0.9999999999999999, 0.0, 0.0);
        FieldVector3D<DerivativeStructure> v2 = createVector(0.0, 0.0, -1.0);

        FieldRotation<DerivativeStructure> quat = new FieldRotation<DerivativeStructure>(u1, u2, v1, v2);
        double q2 = quat.getQ0().getReal() * quat.getQ0().getReal() +
                    quat.getQ1().getReal() * quat.getQ1().getReal() +
                    quat.getQ2().getReal() * quat.getQ2().getReal() +
                    quat.getQ3().getReal() * quat.getQ3().getReal();
        Assert.assertEquals(1.0, q2, 1.0e-14);
        Assert.assertEquals(0.0, FieldVector3D.angle(v1, quat.applyTo(u1)).getReal(), 1.0e-14);
        Assert.assertEquals(0.0, FieldVector3D.angle(v2, quat.applyTo(u2)).getReal(), 1.0e-14);

    }

    @Test
    public void testIssue801() throws MathArithmeticException {
        FieldVector3D<Dfp> u1 = createVector(0.9999988431610581, -0.0015210774290851095, 0.0);
        FieldVector3D<Dfp> u2 = createVector(0.0, 0.0, 1.0);

        FieldVector3D<Dfp> v1 = createVector(0.9999999999999999, 0.0, 0.0);
        FieldVector3D<Dfp> v2 = createVector(0.0, 0.0, -1.0);

        FieldRotation<Dfp> quat = new FieldRotation<Dfp>(u1, u2, v1, v2);
        double q2 = quat.getQ0().getReal() * quat.getQ0().getReal() +
                    quat.getQ1().getReal() * quat.getQ1().getReal() +
                    quat.getQ2().getReal() * quat.getQ2().getReal() +
                    quat.getQ3().getReal() * quat.getQ3().getReal();
        Assert.assertEquals(1.0, q2, 1.0e-14);
        Assert.assertEquals(0.0, FieldVector3D.angle(v1, quat.applyTo(u1)).getReal(), 1.0e-14);
        Assert.assertEquals(0.0, FieldVector3D.angle(v2, quat.applyTo(u2)).getReal(), 1.0e-14);

    }

[193, 'src/test/java', 'org.apache.commons.math3.fraction', 'BigFractionTest', 'testDigitLimitConstructor', 147, 160]

[193, 'src/test/java', 'org.apache.commons.math3.fraction', 'FractionTest', 'testDigitLimitConstructor', 121, 134]

    @Test
    public void testDigitLimitConstructor() throws ConvergenceException {
        assertFraction(2, 5, new BigFraction(0.4, 9));
        assertFraction(2, 5, new BigFraction(0.4, 99));
        assertFraction(2, 5, new BigFraction(0.4, 999));

        assertFraction(3, 5, new BigFraction(0.6152, 9));
        assertFraction(8, 13, new BigFraction(0.6152, 99));
        assertFraction(510, 829, new BigFraction(0.6152, 999));
        assertFraction(769, 1250, new BigFraction(0.6152, 9999));

        // MATH-996
        assertFraction(1, 2, new BigFraction(0.5000000001, 10));
    }

    @Test
    public void testDigitLimitConstructor() throws ConvergenceException  {
        assertFraction(2, 5, new Fraction(0.4,   9));
        assertFraction(2, 5, new Fraction(0.4,  99));
        assertFraction(2, 5, new Fraction(0.4, 999));

        assertFraction(3, 5,      new Fraction(0.6152,    9));
        assertFraction(8, 13,     new Fraction(0.6152,   99));
        assertFraction(510, 829,  new Fraction(0.6152,  999));
        assertFraction(769, 1250, new Fraction(0.6152, 9999));

        // MATH-996
        assertFraction(1, 2, new Fraction(0.5000000001, 10));
    }

[222, 'src/test/java', 'org.apache.commons.math3.linear', 'Array2DRowRealMatrixTest', 'testSerial', 1077, 1081]

[222, 'src/test/java', 'org.apache.commons.math3.linear', 'BlockRealMatrixTest', 'testSerial', 1229, 1233]

    @Test
    public void testSerial()  {
        Array2DRowRealMatrix m = new Array2DRowRealMatrix(testData);
        Assert.assertEquals(m,TestUtils.serializeAndRecover(m));
    }

    @Test
    public void testSerial()  {
        BlockRealMatrix m = new BlockRealMatrix(testData);
        Assert.assertEquals(m,TestUtils.serializeAndRecover(m));
    }

[229, 'src/test/java', 'org.apache.commons.math3.linear', 'Array2DRowRealMatrixTest', 'testFrobeniusNorm', 163, 169]

[229, 'src/test/java', 'org.apache.commons.math3.linear', 'BlockRealMatrixTest', 'testFrobeniusNorm', 164, 170]

    @Test
    public void testFrobeniusNorm() {
        Array2DRowRealMatrix m = new Array2DRowRealMatrix(testData);
        Array2DRowRealMatrix m2 = new Array2DRowRealMatrix(testData2);
        Assert.assertEquals("testData Frobenius norm", FastMath.sqrt(117.0), m.getFrobeniusNorm(), entryTolerance);
        Assert.assertEquals("testData2 Frobenius norm", FastMath.sqrt(52.0), m2.getFrobeniusNorm(), entryTolerance);
    }

    @Test
    public void testFrobeniusNorm() {
        BlockRealMatrix m = new BlockRealMatrix(testData);
        BlockRealMatrix m2 = new BlockRealMatrix(testData2);
        Assert.assertEquals("testData Frobenius norm", FastMath.sqrt(117.0), m.getFrobeniusNorm(), entryTolerance);
        Assert.assertEquals("testData2 Frobenius norm", FastMath.sqrt(52.0), m2.getFrobeniusNorm(), entryTolerance);
    }

[230, 'src/test/java', 'org.apache.commons.math3.stat.descriptive', 'SummaryStatisticsTest', 'testOverrideMeanWithMathClass', 332, 341]

[230, 'src/test/java', 'org.apache.commons.math3.stat.descriptive', 'SummaryStatisticsTest', 'testOverrideGeoMeanWithMathClass', 343, 352]

    @Test
    public void testOverrideMeanWithMathClass() {
        double[] scores = {1, 2, 3, 4};
        SummaryStatistics stats = new SummaryStatistics();
        stats.setMeanImpl(new Mean());
        for(double i : scores) {
          stats.addValue(i);
        }
        Assert.assertEquals((new Mean()).evaluate(scores),stats.getMean(), 0);
    }

    @Test
    public void testOverrideGeoMeanWithMathClass() {
        double[] scores = {1, 2, 3, 4};
        SummaryStatistics stats = new SummaryStatistics();
        stats.setGeoMeanImpl(new GeometricMean());
        for(double i : scores) {
          stats.addValue(i);
        }
        Assert.assertEquals((new GeometricMean()).evaluate(scores),stats.getGeometricMean(), 0);
    }

[255, 'src/test/java', 'org.apache.commons.math3.linear', 'MatrixUtilsTest', 'testCreateRowFieldMatrix', 137, 155]

[255, 'src/test/java', 'org.apache.commons.math3.linear', 'MatrixUtilsTest', 'testCreateColumnFieldMatrix', 175, 194]

    @Test
    public void testCreateRowFieldMatrix() {
        Assert.assertEquals(MatrixUtils.createRowFieldMatrix(asFraction(row)),
                     new Array2DRowFieldMatrix<Fraction>(asFraction(rowMatrix)));
        Assert.assertEquals(MatrixUtils.createRowFieldMatrix(fractionRow),
                     new Array2DRowFieldMatrix<Fraction>(fractionRowMatrix));
        try {
            MatrixUtils.createRowFieldMatrix(new Fraction[] {});  // empty
            Assert.fail("Expecting MathIllegalArgumentException");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }
        try {
            MatrixUtils.createRowFieldMatrix((Fraction[]) null);  // null
            Assert.fail("Expecting NullArgumentException");
        } catch (NullArgumentException ex) {
            // expected
        }
    }

    @Test
    public void testCreateColumnFieldMatrix() {
        Assert.assertEquals(MatrixUtils.createColumnFieldMatrix(asFraction(col)),
                     new Array2DRowFieldMatrix<Fraction>(asFraction(colMatrix)));
        Assert.assertEquals(MatrixUtils.createColumnFieldMatrix(fractionCol),
                     new Array2DRowFieldMatrix<Fraction>(fractionColMatrix));

        try {
            MatrixUtils.createColumnFieldMatrix(new Fraction[] {});  // empty
            Assert.fail("Expecting MathIllegalArgumentException");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }
        try {
            MatrixUtils.createColumnFieldMatrix((Fraction[]) null);  // null
            Assert.fail("Expecting NullArgumentException");
        } catch (NullArgumentException ex) {
            // expected
        }
    }

[274, 'src/test/java', 'org.apache.commons.math3.ode.events', 'ReappearingEventTest', 'testDormandPrince', 36, 42]

[274, 'src/test/java', 'org.apache.commons.math3.ode.events', 'ReappearingEventTest', 'testGragg', 44, 50]

    @Test
    public void testDormandPrince()
        throws DimensionMismatchException, NumberIsTooSmallException,
               MaxCountExceededException, NoBracketingException {
        double tEnd = test(1);
        Assert.assertEquals(10.0, tEnd, 1e-7);
    }

    @Test
    public void testGragg()
        throws DimensionMismatchException, NumberIsTooSmallException,
               MaxCountExceededException, NoBracketingException {
        double tEnd = test(2);
        Assert.assertEquals(10.0, tEnd, 1e-7);
    }

[307, 'src/test/java', 'org.apache.commons.math3.linear', 'RealVectorAbstractTest', 'testSetSubVectorSameType', 395, 408]

[307, 'src/test/java', 'org.apache.commons.math3.linear', 'RealVectorAbstractTest', 'testSetSubVectorMixedType', 410, 423]

    @Test
    public void testSetSubVectorSameType() {
        final double x = getPreferredEntryValue();
        final double[] expected = {x, x, x, 1d, x, 2d, x, x, 3d, x, x, x, 4d, x, x, x};
        final double[] sub = {5d, x, 6d, 7d, 8d};
        final RealVector actual = create(expected);
        final int index = 2;
        actual.setSubVector(index, create(sub));

        for (int i = 0; i < sub.length; i++){
            expected[index + i] = sub[i];
        }
        TestUtils.assertEquals("", expected, actual, 0d);
    }

    @Test
    public void testSetSubVectorMixedType() {
        final double x = getPreferredEntryValue();
        final double[] expected = {x, x, x, 1d, x, 2d, x, x, 3d, x, x, x, 4d, x, x, x};
        final double[] sub = {5d, x, 6d, 7d, 8d};
        final RealVector actual = create(expected);
        final int index = 2;
        actual.setSubVector(index, createAlien(sub));

        for (int i = 0; i < sub.length; i++){
            expected[index + i] = sub[i];
        }
        TestUtils.assertEquals("", expected, actual, 0d);
    }

[336, 'src/test/java', 'org.apache.commons.math3.optimization.direct', 'MultivariateFunctionMappingAdapterTest', 'testStartSimplexInsideRange', 30, 54]

[336, 'src/test/java', 'org.apache.commons.math3.optimization.direct', 'MultivariateFunctionMappingAdapterTest', 'testOptimumOutsideRange', 56, 80]

    @Test
    public void testStartSimplexInsideRange() {

        final BiQuadratic biQuadratic = new BiQuadratic(2.0, 2.5, 1.0, 3.0, 2.0, 3.0);
        final MultivariateFunctionMappingAdapter wrapped =
                new MultivariateFunctionMappingAdapter(biQuadratic,
                                                           biQuadratic.getLower(),
                                                           biQuadratic.getUpper());

        SimplexOptimizer optimizer = new SimplexOptimizer(1e-10, 1e-30);
        optimizer.setSimplex(new NelderMeadSimplex(new double[][] {
            wrapped.boundedToUnbounded(new double[] { 1.5, 2.75 }),
            wrapped.boundedToUnbounded(new double[] { 1.5, 2.95 }),
            wrapped.boundedToUnbounded(new double[] { 1.7, 2.90 })
        }));

        final PointValuePair optimum
            = optimizer.optimize(300, wrapped, GoalType.MINIMIZE,
                                 wrapped.boundedToUnbounded(new double[] { 1.5, 2.25 }));
        final double[] bounded = wrapped.unboundedToBounded(optimum.getPoint());

        Assert.assertEquals(biQuadratic.getBoundedXOptimum(), bounded[0], 2e-7);
        Assert.assertEquals(biQuadratic.getBoundedYOptimum(), bounded[1], 2e-7);

    }

    @Test
    public void testOptimumOutsideRange() {

        final BiQuadratic biQuadratic = new BiQuadratic(4.0, 0.0, 1.0, 3.0, 2.0, 3.0);
        final MultivariateFunctionMappingAdapter wrapped =
                new MultivariateFunctionMappingAdapter(biQuadratic,
                                                           biQuadratic.getLower(),
                                                           biQuadratic.getUpper());

        SimplexOptimizer optimizer = new SimplexOptimizer(1e-10, 1e-30);
        optimizer.setSimplex(new NelderMeadSimplex(new double[][] {
            wrapped.boundedToUnbounded(new double[] { 1.5, 2.75 }),
            wrapped.boundedToUnbounded(new double[] { 1.5, 2.95 }),
            wrapped.boundedToUnbounded(new double[] { 1.7, 2.90 })
        }));

        final PointValuePair optimum
            = optimizer.optimize(100, wrapped, GoalType.MINIMIZE,
                                 wrapped.boundedToUnbounded(new double[] { 1.5, 2.25 }));
        final double[] bounded = wrapped.unboundedToBounded(optimum.getPoint());

        Assert.assertEquals(biQuadratic.getBoundedXOptimum(), bounded[0], 2e-7);
        Assert.assertEquals(biQuadratic.getBoundedYOptimum(), bounded[1], 2e-7);

    }

[354, 'src/test/java', 'org.apache.commons.math3.geometry.euclidean.threed', 'RotationTest', 'testAxisAngleVectorOperator', 83, 114]

[354, 'src/test/java', 'org.apache.commons.math3.geometry.euclidean.threed', 'RotationTest', 'testAxisAngleFrameTransform', 116, 147]

  @Test
  public void testAxisAngleVectorOperator() throws MathIllegalArgumentException {

    Rotation r = new Rotation(new Vector3D(10, 10, 10), 2 * FastMath.PI / 3, RotationConvention.VECTOR_OPERATOR);
    checkVector(r.applyTo(Vector3D.PLUS_I), Vector3D.PLUS_J);
    checkVector(r.applyTo(Vector3D.PLUS_J), Vector3D.PLUS_K);
    checkVector(r.applyTo(Vector3D.PLUS_K), Vector3D.PLUS_I);
    double s = 1 / FastMath.sqrt(3);
    checkVector(r.getAxis(RotationConvention.VECTOR_OPERATOR), new Vector3D( s,  s,  s));
    checkVector(r.getAxis(RotationConvention.FRAME_TRANSFORM), new Vector3D(-s, -s, -s));
    checkAngle(r.getAngle(), 2 * FastMath.PI / 3);

    try {
      new Rotation(new Vector3D(0, 0, 0), 2 * FastMath.PI / 3, RotationConvention.VECTOR_OPERATOR);
      Assert.fail("an exception should have been thrown");
    } catch (MathIllegalArgumentException e) {
    }

    r = new Rotation(Vector3D.PLUS_K, 1.5 * FastMath.PI, RotationConvention.VECTOR_OPERATOR);
    checkVector(r.getAxis(RotationConvention.VECTOR_OPERATOR), new Vector3D(0, 0, -1));
    checkVector(r.getAxis(RotationConvention.FRAME_TRANSFORM), new Vector3D(0, 0, +1));
    checkAngle(r.getAngle(), 0.5 * FastMath.PI);

    r = new Rotation(Vector3D.PLUS_J, FastMath.PI, RotationConvention.VECTOR_OPERATOR);
    checkVector(r.getAxis(RotationConvention.VECTOR_OPERATOR), Vector3D.PLUS_J);
    checkVector(r.getAxis(RotationConvention.FRAME_TRANSFORM), Vector3D.MINUS_J);
    checkAngle(r.getAngle(), FastMath.PI);

    checkVector(Rotation.IDENTITY.getAxis(RotationConvention.VECTOR_OPERATOR), Vector3D.PLUS_I);
    checkVector(Rotation.IDENTITY.getAxis(RotationConvention.FRAME_TRANSFORM), Vector3D.MINUS_I);

  }

  @Test
  public void testAxisAngleFrameTransform() throws MathIllegalArgumentException {

    Rotation r = new Rotation(new Vector3D(10, 10, 10), 2 * FastMath.PI / 3, RotationConvention.FRAME_TRANSFORM);
    checkVector(r.applyTo(Vector3D.PLUS_I), Vector3D.PLUS_K);
    checkVector(r.applyTo(Vector3D.PLUS_J), Vector3D.PLUS_I);
    checkVector(r.applyTo(Vector3D.PLUS_K), Vector3D.PLUS_J);
    double s = 1 / FastMath.sqrt(3);
    checkVector(r.getAxis(RotationConvention.FRAME_TRANSFORM), new Vector3D( s,  s,  s));
    checkVector(r.getAxis(RotationConvention.VECTOR_OPERATOR), new Vector3D(-s, -s, -s));
    checkAngle(r.getAngle(), 2 * FastMath.PI / 3);

    try {
      new Rotation(new Vector3D(0, 0, 0), 2 * FastMath.PI / 3, RotationConvention.FRAME_TRANSFORM);
      Assert.fail("an exception should have been thrown");
    } catch (MathIllegalArgumentException e) {
    }

    r = new Rotation(Vector3D.PLUS_K, 1.5 * FastMath.PI, RotationConvention.FRAME_TRANSFORM);
    checkVector(r.getAxis(RotationConvention.FRAME_TRANSFORM), new Vector3D(0, 0, -1));
    checkVector(r.getAxis(RotationConvention.VECTOR_OPERATOR), new Vector3D(0, 0, +1));
    checkAngle(r.getAngle(), 0.5 * FastMath.PI);

    r = new Rotation(Vector3D.PLUS_J, FastMath.PI, RotationConvention.FRAME_TRANSFORM);
    checkVector(r.getAxis(RotationConvention.FRAME_TRANSFORM), Vector3D.PLUS_J);
    checkVector(r.getAxis(RotationConvention.VECTOR_OPERATOR), Vector3D.MINUS_J);
    checkAngle(r.getAngle(), FastMath.PI);

    checkVector(Rotation.IDENTITY.getAxis(RotationConvention.FRAME_TRANSFORM), Vector3D.MINUS_I);
    checkVector(Rotation.IDENTITY.getAxis(RotationConvention.VECTOR_OPERATOR), Vector3D.PLUS_I);

  }

[392, 'src/test/java', 'org.apache.commons.math3.linear', 'QRSolverTest', 'testOverdetermined', 145, 167]

[392, 'src/test/java', 'org.apache.commons.math3.linear', 'RRQRSolverTest', 'testOverdetermined', 148, 170]

    @Test
    public void testOverdetermined() {
        final Random r    = new Random(5559252868205245l);
        int          p    = (7 * BlockRealMatrix.BLOCK_SIZE) / 4;
        int          q    = (5 * BlockRealMatrix.BLOCK_SIZE) / 4;
        RealMatrix   a    = createTestMatrix(r, p, q);
        RealMatrix   xRef = createTestMatrix(r, q, BlockRealMatrix.BLOCK_SIZE + 3);

        // build a perturbed system: A.X + noise = B
        RealMatrix b = a.multiply(xRef);
        final double noise = 0.001;
        b.walkInOptimizedOrder(new DefaultRealMatrixChangingVisitor() {
            @Override
            public double visit(int row, int column, double value) {
                return value * (1.0 + noise * (2 * r.nextDouble() - 1));
            }
        });

        // despite perturbation, the least square solution should be pretty good
        RealMatrix x = new QRDecomposition(a).getSolver().solve(b);
        Assert.assertEquals(0, x.subtract(xRef).getNorm(), 0.01 * noise * p * q);

    }

    @Test
    public void testOverdetermined() {
        final Random r    = new Random(5559252868205245l);
        int          p    = (7 * BlockRealMatrix.BLOCK_SIZE) / 4;
        int          q    = (5 * BlockRealMatrix.BLOCK_SIZE) / 4;
        RealMatrix   a    = createTestMatrix(r, p, q);
        RealMatrix   xRef = createTestMatrix(r, q, BlockRealMatrix.BLOCK_SIZE + 3);

        // build a perturbed system: A.X + noise = B
        RealMatrix b = a.multiply(xRef);
        final double noise = 0.001;
        b.walkInOptimizedOrder(new DefaultRealMatrixChangingVisitor() {
            @Override
            public double visit(int row, int column, double value) {
                return value * (1.0 + noise * (2 * r.nextDouble() - 1));
            }
        });

        // despite perturbation, the least square solution should be pretty good
        RealMatrix x = new RRQRDecomposition(a).getSolver().solve(b);
        Assert.assertEquals(0, x.subtract(xRef).getNorm(), 0.01 * noise * p * q);

    }

[401, 'src/test/java', 'org.apache.commons.math3.linear', 'HessenbergTransformerTest', 'testRandomDataNormalDistribution', 107, 127]

[401, 'src/test/java', 'org.apache.commons.math3.linear', 'SchurTransformerTest', 'testRandomDataNormalDistribution', 111, 131]

    @Test
    public void testRandomDataNormalDistribution() {
        for (int run = 0; run < 100; run++) {
            Random r = new Random(System.currentTimeMillis());
            NormalDistribution dist = new NormalDistribution(0.0, r.nextDouble() * 5);

            // matrix size
            int size = r.nextInt(20) + 4;

            double[][] data = new double[size][size];
            for (int i = 0; i < size; i++) {
                for (int j = 0; j < size; j++) {
                    data[i][j] = dist.sample();
                }
            }

            RealMatrix m = MatrixUtils.createRealMatrix(data);
            RealMatrix h = checkAEqualPHPt(m);
            checkHessenbergForm(h);
        }
    }

    @Test
    public void testRandomDataNormalDistribution() {
        for (int run = 0; run < 100; run++) {
            Random r = new Random(System.currentTimeMillis());
            NormalDistribution dist = new NormalDistribution(0.0, r.nextDouble() * 5);

            // matrix size
            int size = r.nextInt(20) + 4;

            double[][] data = new double[size][size];
            for (int i = 0; i < size; i++) {
                for (int j = 0; j < size; j++) {
                    data[i][j] = dist.sample();
                }
            }

            RealMatrix m = MatrixUtils.createRealMatrix(data);
            RealMatrix s = checkAEqualPTPt(m);
            checkSchurForm(s);
        }
    }

[423, 'src/test/java', 'org.apache.commons.math3.linear', 'Array2DRowRealMatrixTest', 'testGetEntry', 405, 415]

[423, 'src/test/java', 'org.apache.commons.math3.linear', 'BlockRealMatrixTest', 'testGetEntry', 450, 460]

    @Test
    public void testGetEntry() {
        RealMatrix m = new Array2DRowRealMatrix(testData);
        Assert.assertEquals("get entry",m.getEntry(0,1),2d,entryTolerance);
        try {
            m.getEntry(10, 4);
            Assert.fail("Expecting OutOfRangeException");
        } catch (OutOfRangeException ex) {
            // expected
        }
    }

    @Test
    public void testGetEntry() {
        RealMatrix m = new BlockRealMatrix(testData);
        Assert.assertEquals("get entry",m.getEntry(0,1),2d,entryTolerance);
        try {
            m.getEntry(10, 4);
            Assert.fail ("Expecting OutOfRangeException");
        } catch (OutOfRangeException ex) {
            // expected
        }
    }

[424, 'src/test/java', 'org.apache.commons.math3.fraction', 'BigFractionFormatTest', 'testFormatImproper', 83, 92]

[424, 'src/test/java', 'org.apache.commons.math3.fraction', 'FractionFormatTest', 'testFormatImproper', 81, 90]

    @Test
    public void testFormatImproper() {
        BigFraction c = new BigFraction(5, 3);

        String actual = properFormat.format(c);
        Assert.assertEquals("1 2 / 3", actual);

        actual = improperFormat.format(c);
        Assert.assertEquals("5 / 3", actual);
    }

    @Test
    public void testFormatImproper() {
        Fraction c = new Fraction(5, 3);

        String actual = properFormat.format(c);
        Assert.assertEquals("1 2 / 3", actual);

        actual = improperFormat.format(c);
        Assert.assertEquals("5 / 3", actual);
    }

[425, 'src/test/java', 'org.apache.commons.math3.fraction', 'BigFractionFormatTest', 'testWholeFormat', 311, 321]

[425, 'src/test/java', 'org.apache.commons.math3.fraction', 'FractionFormatTest', 'testWholeFormat', 330, 340]

    @Test
    public void testWholeFormat() {
        ProperBigFractionFormat format = (ProperBigFractionFormat)properFormat;

        NumberFormat old = format.getWholeFormat();
        NumberFormat nf = NumberFormat.getInstance();
        nf.setParseIntegerOnly(true);
        format.setWholeFormat(nf);
        Assert.assertEquals(nf, format.getWholeFormat());
        format.setWholeFormat(old);
    }

    @Test
    public void testWholeFormat() {
        ProperFractionFormat format = (ProperFractionFormat)properFormat;

        NumberFormat old = format.getWholeFormat();
        NumberFormat nf = NumberFormat.getInstance();
        nf.setParseIntegerOnly(true);
        format.setWholeFormat(nf);
        Assert.assertEquals(nf, format.getWholeFormat());
        format.setWholeFormat(old);
    }

[442, 'src/test/java', 'org.apache.commons.math3.analysis.polynomials', 'PolynomialsUtilsTest', 'testFirstChebyshevPolynomials', 33, 45]

[442, 'src/test/java', 'org.apache.commons.math3.analysis.polynomials', 'PolynomialsUtilsTest', 'testFirstHermitePolynomials', 94, 106]

    @Test
    public void testFirstChebyshevPolynomials() {
        checkPolynomial(PolynomialsUtils.createChebyshevPolynomial(3), "-3 x + 4 x^3");
        checkPolynomial(PolynomialsUtils.createChebyshevPolynomial(2), "-1 + 2 x^2");
        checkPolynomial(PolynomialsUtils.createChebyshevPolynomial(1), "x");
        checkPolynomial(PolynomialsUtils.createChebyshevPolynomial(0), "1");

        checkPolynomial(PolynomialsUtils.createChebyshevPolynomial(7), "-7 x + 56 x^3 - 112 x^5 + 64 x^7");
        checkPolynomial(PolynomialsUtils.createChebyshevPolynomial(6), "-1 + 18 x^2 - 48 x^4 + 32 x^6");
        checkPolynomial(PolynomialsUtils.createChebyshevPolynomial(5), "5 x - 20 x^3 + 16 x^5");
        checkPolynomial(PolynomialsUtils.createChebyshevPolynomial(4), "1 - 8 x^2 + 8 x^4");

    }

    @Test
    public void testFirstHermitePolynomials() {
        checkPolynomial(PolynomialsUtils.createHermitePolynomial(3), "-12 x + 8 x^3");
        checkPolynomial(PolynomialsUtils.createHermitePolynomial(2), "-2 + 4 x^2");
        checkPolynomial(PolynomialsUtils.createHermitePolynomial(1), "2 x");
        checkPolynomial(PolynomialsUtils.createHermitePolynomial(0), "1");

        checkPolynomial(PolynomialsUtils.createHermitePolynomial(7), "-1680 x + 3360 x^3 - 1344 x^5 + 128 x^7");
        checkPolynomial(PolynomialsUtils.createHermitePolynomial(6), "-120 + 720 x^2 - 480 x^4 + 64 x^6");
        checkPolynomial(PolynomialsUtils.createHermitePolynomial(5), "120 x - 160 x^3 + 32 x^5");
        checkPolynomial(PolynomialsUtils.createHermitePolynomial(4), "12 - 48 x^2 + 16 x^4");

    }

[455, 'src/test/java', 'org.apache.commons.math3.distribution', 'TriangularDistributionTest', 'testGetLowerBound', 139, 143]

[455, 'src/test/java', 'org.apache.commons.math3.distribution', 'UniformRealDistributionTest', 'testGetLowerBound', 71, 75]

    @Test
    public void testGetLowerBound() {
        TriangularDistribution distribution = makeDistribution();
        Assert.assertEquals(-3.0, distribution.getSupportLowerBound(), 0);
    }

    @Test
    public void testGetLowerBound() {
        UniformRealDistribution distribution = makeDistribution();
        Assert.assertEquals(-0.5, distribution.getSupportLowerBound(), 0);
    }

[463, 'src/test/java', 'org.apache.commons.math3.fitting', 'GaussianCurveFitterTest', 'testFit01', 181, 189]

[463, 'src/test/java', 'org.apache.commons.math3.fitting', 'GaussianCurveFitterTest', 'testFit07', 294, 302]

    @Test
    public void testFit01() {
        GaussianCurveFitter fitter = GaussianCurveFitter.create();
        double[] parameters = fitter.fit(createDataset(DATASET1).toList());

        Assert.assertEquals(3496978.1837704973, parameters[0], 1e-4);
        Assert.assertEquals(4.054933085999146, parameters[1], 1e-4);
        Assert.assertEquals(0.015039355620304326, parameters[2], 1e-4);
    }

    @Test
    public void testFit07() {
        GaussianCurveFitter fitter = GaussianCurveFitter.create();
        double[] parameters = fitter.fit(createDataset(DATASET5).toList());

        Assert.assertEquals(3514384.729342235, parameters[0], 1e-4);
        Assert.assertEquals(4.054970307455625, parameters[1], 1e-4);
        Assert.assertEquals(0.015029412832160017, parameters[2], 1e-4);
    }

[479, 'src/test/java', 'org.apache.commons.math3.stat.regression', 'MultipleLinearRegressionAbstractTest', 'testNewSampleInvalidData', 113, 117]

[479, 'src/test/java', 'org.apache.commons.math3.stat.regression', 'MultipleLinearRegressionAbstractTest', 'testNewSampleInsufficientData', 119, 123]

    @Test(expected=IllegalArgumentException.class)
    public void testNewSampleInvalidData() {
        double[] data = new double[] {1, 2, 3, 4};
        createRegression().newSampleData(data, 2, 3);
    }

    @Test(expected=IllegalArgumentException.class)
    public void testNewSampleInsufficientData() {
        double[] data = new double[] {1, 2, 3, 4};
        createRegression().newSampleData(data, 1, 3);
    }

[494, 'src/test/java', 'org.apache.commons.math3.fraction', 'BigFractionTest', 'testDoubleValue', 208, 215]

[494, 'src/test/java', 'org.apache.commons.math3.fraction', 'FractionTest', 'testDoubleValue', 188, 195]

    @Test
    public void testDoubleValue() {
        BigFraction first = new BigFraction(1, 2);
        BigFraction second = new BigFraction(1, 3);

        Assert.assertEquals(0.5, first.doubleValue(), 0.0);
        Assert.assertEquals(1.0 / 3.0, second.doubleValue(), 0.0);
    }

    @Test
    public void testDoubleValue() {
        Fraction first = new Fraction(1, 2);
        Fraction second = new Fraction(1, 3);

        Assert.assertEquals(0.5, first.doubleValue(), 0.0);
        Assert.assertEquals(1.0 / 3.0, second.doubleValue(), 0.0);
    }

[506, 'src/test/java', 'org.apache.commons.math3.analysis.differentiation', 'DSCompilerTest', 'testIncompatibleParams', 132, 135]

[506, 'src/test/java', 'org.apache.commons.math3.analysis.differentiation', 'DSCompilerTest', 'testIncompatibleOrder', 137, 140]

    @Test(expected=DimensionMismatchException.class)
    public void testIncompatibleParams() {
        DSCompiler.getCompiler(3, 2).checkCompatibility(DSCompiler.getCompiler(4, 2));
    }

    @Test(expected=DimensionMismatchException.class)
    public void testIncompatibleOrder() {
        DSCompiler.getCompiler(3, 3).checkCompatibility(DSCompiler.getCompiler(3, 2));
    }

[535, 'src/test/java', 'org.apache.commons.math3.analysis.solvers', 'NewtonRaphsonSolverTest', 'testSinZero', 33, 46]

[535, 'src/test/java', 'org.apache.commons.math3.analysis.solvers', 'NewtonSolverTest', 'testSinZero', 38, 51]

    @Test
    public void testSinZero() {
        UnivariateDifferentiableFunction f = new Sin();
        double result;

        NewtonRaphsonSolver solver = new NewtonRaphsonSolver();
        result = solver.solve(100, f, 3, 4);
        Assert.assertEquals(result, FastMath.PI, solver.getAbsoluteAccuracy());

        result = solver.solve(100, f, 1, 4);
        Assert.assertEquals(result, FastMath.PI, solver.getAbsoluteAccuracy());

        Assert.assertTrue(solver.getEvaluations() > 0);
    }

    @Test
    public void testSinZero() {
        DifferentiableUnivariateFunction f = new Sin();
        double result;

        NewtonSolver solver = new NewtonSolver();
        result = solver.solve(100, f, 3, 4);
        Assert.assertEquals(result, FastMath.PI, solver.getAbsoluteAccuracy());

        result = solver.solve(100, f, 1, 4);
        Assert.assertEquals(result, FastMath.PI, solver.getAbsoluteAccuracy());

        Assert.assertTrue(solver.getEvaluations() > 0);
    }

[538, 'src/test/java', 'org.apache.commons.math3.geometry.euclidean.threed', 'FieldRotationDSTest', 'testApplyInverseTo', 1133, 1185]

[538, 'src/test/java', 'org.apache.commons.math3.geometry.euclidean.threed', 'FieldRotationDfpTest', 'testApplyInverseTo', 900, 950]

    @Test
    public void testApplyInverseTo() throws MathIllegalArgumentException {

        DerivativeStructure[] in      = new DerivativeStructure[3];
        DerivativeStructure[] out     = new DerivativeStructure[3];
        DerivativeStructure[] rebuilt = new DerivativeStructure[3];
        FieldRotation<DerivativeStructure> r = new FieldRotation<DerivativeStructure>(createVector(2, -3, 5),
                                                                                      createAngle(1.7),
                                                                                      RotationConvention.VECTOR_OPERATOR);
        for (double lambda = 0; lambda < 6.2; lambda += 0.2) {
            for (double phi = -1.55; phi < 1.55; phi += 0.2) {
                FieldVector3D<DerivativeStructure> u = createVector(FastMath.cos(lambda) * FastMath.cos(phi),
                                          FastMath.sin(lambda) * FastMath.cos(phi),
                                          FastMath.sin(phi));
                r.applyInverseTo(r.applyTo(u));
                checkVector(u, r.applyInverseTo(r.applyTo(u)));
                checkVector(u, r.applyTo(r.applyInverseTo(u)));
                in[0] = u.getX();
                in[1] = u.getY();
                in[2] = u.getZ();
                r.applyTo(in, out);
                r.applyInverseTo(out, rebuilt);
                Assert.assertEquals(in[0].getReal(), rebuilt[0].getReal(), 1.0e-12);
                Assert.assertEquals(in[1].getReal(), rebuilt[1].getReal(), 1.0e-12);
                Assert.assertEquals(in[2].getReal(), rebuilt[2].getReal(), 1.0e-12);
            }
        }

        r = createRotation(1, 0, 0, 0, false);
        for (double lambda = 0; lambda < 6.2; lambda += 0.2) {
            for (double phi = -1.55; phi < 1.55; phi += 0.2) {
                FieldVector3D<DerivativeStructure> u = createVector(FastMath.cos(lambda) * FastMath.cos(phi),
                                          FastMath.sin(lambda) * FastMath.cos(phi),
                                          FastMath.sin(phi));
                checkVector(u, r.applyInverseTo(r.applyTo(u)));
                checkVector(u, r.applyTo(r.applyInverseTo(u)));
            }
        }

        r = new FieldRotation<DerivativeStructure>(createVector(0, 0, 1),
                                                   createAngle(FastMath.PI),
                                                   RotationConvention.VECTOR_OPERATOR);
        for (double lambda = 0; lambda < 6.2; lambda += 0.2) {
            for (double phi = -1.55; phi < 1.55; phi += 0.2) {
                FieldVector3D<DerivativeStructure> u = createVector(FastMath.cos(lambda) * FastMath.cos(phi),
                                          FastMath.sin(lambda) * FastMath.cos(phi),
                                          FastMath.sin(phi));
                checkVector(u, r.applyInverseTo(r.applyTo(u)));
                checkVector(u, r.applyTo(r.applyInverseTo(u)));
            }
        }

    }

    @Test
    public void testApplyInverseTo() throws MathIllegalArgumentException {

        Dfp[] in      = new Dfp[3];
        Dfp[] out     = new Dfp[3];
        Dfp[] rebuilt = new Dfp[3];
        FieldRotation<Dfp> r = new FieldRotation<Dfp>(createVector(2, -3, 5),
                                                      createAngle(1.7),
                                                      RotationConvention.VECTOR_OPERATOR);
        for (double lambda = 0; lambda < 6.2; lambda += 0.2) {
            for (double phi = -1.55; phi < 1.55; phi += 0.2) {
                FieldVector3D<Dfp> u = createVector(FastMath.cos(lambda) * FastMath.cos(phi),
                                          FastMath.sin(lambda) * FastMath.cos(phi),
                                          FastMath.sin(phi));
                r.applyInverseTo(r.applyTo(u));
                checkVector(u, r.applyInverseTo(r.applyTo(u)));
                checkVector(u, r.applyTo(r.applyInverseTo(u)));
                in[0] = u.getX();
                in[1] = u.getY();
                in[2] = u.getZ();
                r.applyTo(in, out);
                r.applyInverseTo(out, rebuilt);
                Assert.assertEquals(in[0].getReal(), rebuilt[0].getReal(), 1.0e-12);
                Assert.assertEquals(in[1].getReal(), rebuilt[1].getReal(), 1.0e-12);
                Assert.assertEquals(in[2].getReal(), rebuilt[2].getReal(), 1.0e-12);
            }
        }

        r = createRotation(1, 0, 0, 0, false);
        for (double lambda = 0; lambda < 6.2; lambda += 0.2) {
            for (double phi = -1.55; phi < 1.55; phi += 0.2) {
                FieldVector3D<Dfp> u = createVector(FastMath.cos(lambda) * FastMath.cos(phi),
                                          FastMath.sin(lambda) * FastMath.cos(phi),
                                          FastMath.sin(phi));
                checkVector(u, r.applyInverseTo(r.applyTo(u)));
                checkVector(u, r.applyTo(r.applyInverseTo(u)));
            }
        }

        r = new FieldRotation<Dfp>(createVector(0, 0, 1), createAngle(FastMath.PI), RotationConvention.VECTOR_OPERATOR);
        for (double lambda = 0; lambda < 6.2; lambda += 0.2) {
            for (double phi = -1.55; phi < 1.55; phi += 0.2) {
                FieldVector3D<Dfp> u = createVector(FastMath.cos(lambda) * FastMath.cos(phi),
                                          FastMath.sin(lambda) * FastMath.cos(phi),
                                          FastMath.sin(phi));
                checkVector(u, r.applyInverseTo(r.applyTo(u)));
                checkVector(u, r.applyTo(r.applyInverseTo(u)));
            }
        }

    }

[565, 'src/test/java', 'org.apache.commons.math3.complex', 'ComplexTest', 'testFloatingPointEqualsWithAllowedDeltaNaN', 567, 574]

[565, 'src/test/java', 'org.apache.commons.math3.complex', 'ComplexTest', 'testFloatingPointEqualsWithRelativeToleranceNaN', 588, 595]

    @Test
    public void testFloatingPointEqualsWithAllowedDeltaNaN() {
        final Complex x = new Complex(0, Double.NaN);
        final Complex y = new Complex(Double.NaN, 0);
        Assert.assertFalse(Complex.equals(x, Complex.ZERO, 0.1));
        Assert.assertFalse(Complex.equals(x, x, 0.1));
        Assert.assertFalse(Complex.equals(x, y, 0.1));
    }

    @Test
    public void testFloatingPointEqualsWithRelativeToleranceNaN() {
        final Complex x = new Complex(0, Double.NaN);
        final Complex y = new Complex(Double.NaN, 0);
        Assert.assertFalse(Complex.equalsWithRelativeTolerance(x, Complex.ZERO, 0.1));
        Assert.assertFalse(Complex.equalsWithRelativeTolerance(x, x, 0.1));
        Assert.assertFalse(Complex.equalsWithRelativeTolerance(x, y, 0.1));
    }

[574, 'src/test/java', 'org.apache.commons.math3.linear', 'BlockFieldMatrixTest', 'testAddFail', 199, 209]

[574, 'src/test/java', 'org.apache.commons.math3.linear', 'FieldMatrixImplTest', 'testAddFail', 140, 150]

    @Test
    public void testAddFail() {
        BlockFieldMatrix<Fraction> m = new BlockFieldMatrix<Fraction>(testData);
        BlockFieldMatrix<Fraction> m2 = new BlockFieldMatrix<Fraction>(testData2);
        try {
            m.add(m2);
            Assert.fail("MathIllegalArgumentException expected");
        } catch (MathIllegalArgumentException ex) {
            // ignored
        }
    }

    @Test
    public void testAddFail() {
        Array2DRowFieldMatrix<Fraction> m = new Array2DRowFieldMatrix<Fraction>(testData);
        Array2DRowFieldMatrix<Fraction> m2 = new Array2DRowFieldMatrix<Fraction>(testData2);
        try {
            m.add(m2);
            Assert.fail("MathIllegalArgumentException expected");
        } catch (MathIllegalArgumentException ex) {
            // ignored
        }
    }

[588, 'src/test/java', 'org.apache.commons.math3.util', 'ArithmeticUtilsTest', 'testPowNegativeLong', 534, 542]

[588, 'src/test/java', 'org.apache.commons.math3.util', 'MathUtilsTest', 'testSignLong', 243, 248]

    @Test
    public void testPowNegativeLong() {
        final long base = -21;

        Assert.assertEquals(-154472377739119461L,
                            ArithmeticUtils.pow(base, 13));
        Assert.assertEquals(3243919932521508681L,
                            ArithmeticUtils.pow(base, 14));
    }

    @Test
    public void testSignLong() {
        final long one = 1L;
        Assert.assertEquals(1L, MathUtils.copySign(one, 2L));
        Assert.assertEquals(-1L, MathUtils.copySign(one, -2L));
    }

[591, 'src/test/java', 'org.apache.commons.math3.geometry.euclidean.threed', 'Vector3DFormatAbstractTest', 'testSimpleNoDecimals', 46, 52]

[591, 'src/test/java', 'org.apache.commons.math3.geometry.euclidean.threed', 'Vector3DFormatAbstractTest', 'testNonDefaultSetting', 114, 120]

    @Test
    public void testSimpleNoDecimals() {
        Vector3D c = new Vector3D(1, 1, 1);
        String expected = "{1; 1; 1}";
        String actual = vector3DFormat.format(c);
        Assert.assertEquals(expected, actual);
    }

    @Test
    public void testNonDefaultSetting() {
        Vector3D c = new Vector3D(1, 1, 1);
        String expected = "[1 : 1 : 1]";
        String actual = vector3DFormatSquare.format(c);
        Assert.assertEquals(expected, actual);
    }

[595, 'src/test/java', 'org.apache.commons.math3.analysis.interpolation', 'DividedDifferenceInterpolatorTest', 'testExpm1Function', 80, 109]

[595, 'src/test/java', 'org.apache.commons.math3.analysis.interpolation', 'NevilleInterpolatorTest', 'testExpm1Function', 80, 109]

    @Test
    public void testExpm1Function() {
        UnivariateFunction f = new Expm1();
        UnivariateInterpolator interpolator = new DividedDifferenceInterpolator();
        double x[], y[], z, expected, result, tolerance;

        // 5 interpolating points on interval [-1, 1]
        int n = 5;
        double min = -1.0, max = 1.0;
        x = new double[n];
        y = new double[n];
        for (int i = 0; i < n; i++) {
            x[i] = min + i * (max - min) / n;
            y[i] = f.value(x[i]);
        }
        double derivativebound = FastMath.E;
        UnivariateFunction p = interpolator.interpolate(x, y);

        z = 0.0; expected = f.value(z); result = p.value(z);
        tolerance = FastMath.abs(derivativebound * partialerror(x, z));
        Assert.assertEquals(expected, result, tolerance);

        z = 0.5; expected = f.value(z); result = p.value(z);
        tolerance = FastMath.abs(derivativebound * partialerror(x, z));
        Assert.assertEquals(expected, result, tolerance);

        z = -0.5; expected = f.value(z); result = p.value(z);
        tolerance = FastMath.abs(derivativebound * partialerror(x, z));
        Assert.assertEquals(expected, result, tolerance);
    }

    @Test
    public void testExpm1Function() {
        UnivariateFunction f = new Expm1();
        UnivariateInterpolator interpolator = new NevilleInterpolator();
        double x[], y[], z, expected, result, tolerance;

        // 5 interpolating points on interval [-1, 1]
        int n = 5;
        double min = -1.0, max = 1.0;
        x = new double[n];
        y = new double[n];
        for (int i = 0; i < n; i++) {
            x[i] = min + i * (max - min) / n;
            y[i] = f.value(x[i]);
        }
        double derivativebound = FastMath.E;
        UnivariateFunction p = interpolator.interpolate(x, y);

        z = 0.0; expected = f.value(z); result = p.value(z);
        tolerance = FastMath.abs(derivativebound * partialerror(x, z));
        Assert.assertEquals(expected, result, tolerance);

        z = 0.5; expected = f.value(z); result = p.value(z);
        tolerance = FastMath.abs(derivativebound * partialerror(x, z));
        Assert.assertEquals(expected, result, tolerance);

        z = -0.5; expected = f.value(z); result = p.value(z);
        tolerance = FastMath.abs(derivativebound * partialerror(x, z));
        Assert.assertEquals(expected, result, tolerance);
    }

[601, 'src/test/java', 'org.apache.commons.math3.linear', 'SparseFieldVectorTest', 'testWalkInDefaultOrderPreservingVisitor2', 287, 335]

[601, 'src/test/java', 'org.apache.commons.math3.linear', 'SparseFieldVectorTest', 'testWalkInOptimizedOrderPreservingVisitor2', 414, 462]

    @Test
    public void testWalkInDefaultOrderPreservingVisitor2() {
        final SparseFieldVector<Fraction> v = create(5);
        final FieldVectorPreservingVisitor<Fraction> visitor;
        visitor = new FieldVectorPreservingVisitor<Fraction>() {

            public void visit(int index, Fraction value) {
                // Do nothing
            }

            public void start(int dimension, int start, int end) {
                // Do nothing
            }

            public Fraction end() {
                return Fraction.ZERO;
            }
        };
        try {
            v.walkInDefaultOrder(visitor, -1, 4);
            Assert.fail();
        } catch (OutOfRangeException e) {
            // Expected behavior
        }
        try {
            v.walkInDefaultOrder(visitor, 5, 4);
            Assert.fail();
        } catch (OutOfRangeException e) {
            // Expected behavior
        }
        try {
            v.walkInDefaultOrder(visitor, 0, -1);
            Assert.fail();
        } catch (OutOfRangeException e) {
            // Expected behavior
        }
        try {
            v.walkInDefaultOrder(visitor, 0, 5);
            Assert.fail();
        } catch (OutOfRangeException e) {
            // Expected behavior
        }
        try {
            v.walkInDefaultOrder(visitor, 4, 0);
            Assert.fail();
        } catch (NumberIsTooSmallException e) {
            // Expected behavior
        }
    }

    @Test
    public void testWalkInOptimizedOrderPreservingVisitor2() {
        final SparseFieldVector<Fraction> v = create(5);
        final FieldVectorPreservingVisitor<Fraction> visitor;
        visitor = new FieldVectorPreservingVisitor<Fraction>() {

            public void visit(int index, Fraction value) {
                // Do nothing
            }

            public void start(int dimension, int start, int end) {
                // Do nothing
            }

            public Fraction end() {
                return Fraction.ZERO;
            }
        };
        try {
            v.walkInOptimizedOrder(visitor, -1, 4);
            Assert.fail();
        } catch (OutOfRangeException e) {
            // Expected behavior
        }
        try {
            v.walkInOptimizedOrder(visitor, 5, 4);
            Assert.fail();
        } catch (OutOfRangeException e) {
            // Expected behavior
        }
        try {
            v.walkInOptimizedOrder(visitor, 0, -1);
            Assert.fail();
        } catch (OutOfRangeException e) {
            // Expected behavior
        }
        try {
            v.walkInOptimizedOrder(visitor, 0, 5);
            Assert.fail();
        } catch (OutOfRangeException e) {
            // Expected behavior
        }
        try {
            v.walkInOptimizedOrder(visitor, 4, 0);
            Assert.fail();
        } catch (NumberIsTooSmallException e) {
            // Expected behavior
        }
    }

[610, 'src/test/java', 'org.apache.commons.math3.analysis.function', 'SincTest', 'testDerivatives1Dot2Unnormalized', 82, 91]

[610, 'src/test/java', 'org.apache.commons.math3.analysis.function', 'SincTest', 'testDerivatives1Dot2Normalized', 93, 102]

   @Test
   public void testDerivatives1Dot2Unnormalized() {
       DerivativeStructure s = new Sinc(false).value(new DerivativeStructure(1, 5, 0, 1.2));
       Assert.assertEquals( 0.77669923830602195806, s.getPartialDerivative(0), 1.0e-16);
       Assert.assertEquals(-0.34528456985779031701, s.getPartialDerivative(1), 1.0e-16);
       Assert.assertEquals(-0.2012249552097047631,  s.getPartialDerivative(2), 1.0e-16);
       Assert.assertEquals( 0.2010975926270339262,  s.getPartialDerivative(3), 4.0e-16);
       Assert.assertEquals( 0.106373929549242204,   s.getPartialDerivative(4), 1.0e-15);
       Assert.assertEquals(-0.1412599110579478695,  s.getPartialDerivative(5), 3.0e-15);
   }

   @Test
   public void testDerivatives1Dot2Normalized() {
       DerivativeStructure s = new Sinc(true).value(new DerivativeStructure(1, 5, 0, 1.2));
       Assert.assertEquals(-0.15591488063143983888, s.getPartialDerivative(0), 6.0e-17);
       Assert.assertEquals(-0.54425176145292298767, s.getPartialDerivative(1), 2.0e-16);
       Assert.assertEquals(2.4459044611635856107,   s.getPartialDerivative(2), 9.0e-16);
       Assert.assertEquals(0.5391369206235909586,   s.getPartialDerivative(3), 7.0e-16);
       Assert.assertEquals(-16.984649869728849865,  s.getPartialDerivative(4), 8.0e-15);
       Assert.assertEquals(5.0980327462666316586,   s.getPartialDerivative(5), 9.0e-15);
   }

[616, 'src/test/java', 'org.apache.commons.math3.dfp', 'DfpTest', 'testIntConstructor', 94, 103]

[616, 'src/test/java', 'org.apache.commons.math3.dfp', 'DfpTest', 'testLongConstructor', 105, 114]

    @Test
    public void testIntConstructor() {
        Assert.assertEquals("0.", new Dfp(field, 0).toString());
        Assert.assertEquals("1.", new Dfp(field, 1).toString());
        Assert.assertEquals("-1.", new Dfp(field, -1).toString());
        Assert.assertEquals("1234567890.", new Dfp(field, 1234567890).toString());
        Assert.assertEquals("-1234567890.", new Dfp(field, -1234567890).toString());
        Assert.assertEquals("-2147483648.", new Dfp(field, Integer.MIN_VALUE).toString());
        Assert.assertEquals("2147483647.", new Dfp(field, Integer.MAX_VALUE).toString());
    }

    @Test
    public void testLongConstructor() {
        Assert.assertEquals("0.", new Dfp(field, 0l).toString());
        Assert.assertEquals("1.", new Dfp(field, 1l).toString());
        Assert.assertEquals("-1.", new Dfp(field, -1l).toString());
        Assert.assertEquals("1234567890.", new Dfp(field, 1234567890l).toString());
        Assert.assertEquals("-1234567890.", new Dfp(field, -1234567890l).toString());
        Assert.assertEquals("-9223372036854775808.", new Dfp(field, Long.MIN_VALUE).toString());
        Assert.assertEquals("9223372036854775807.", new Dfp(field, Long.MAX_VALUE).toString());
    }

[617, 'src/test/java', 'org.apache.commons.math3.geometry.euclidean.threed', 'FieldRotationDSTest', 'testVectorOnePair', 386, 403]

[617, 'src/test/java', 'org.apache.commons.math3.geometry.euclidean.threed', 'FieldRotationDfpTest', 'testVectorOnePair', 226, 243]

    @Test
    public void testVectorOnePair() throws MathArithmeticException {

        FieldVector3D<DerivativeStructure> u = createVector(3, 2, 1);
        FieldVector3D<DerivativeStructure> v = createVector(-4, 2, 2);
        FieldRotation<DerivativeStructure> r = new FieldRotation<DerivativeStructure>(u, v);
        checkVector(r.applyTo(u.scalarMultiply(v.getNorm())), v.scalarMultiply(u.getNorm()));

        checkAngle(new FieldRotation<DerivativeStructure>(u, u.negate()).getAngle(), FastMath.PI);

        try {
            new FieldRotation<DerivativeStructure>(u, createVector(0, 0, 0));
            Assert.fail("an exception should have been thrown");
        } catch (MathArithmeticException e) {
            // expected behavior
        }

    }

    @Test
    public void testVectorOnePair() throws MathArithmeticException {

        FieldVector3D<Dfp> u = createVector(3, 2, 1);
        FieldVector3D<Dfp> v = createVector(-4, 2, 2);
        FieldRotation<Dfp> r = new FieldRotation<Dfp>(u, v);
        checkVector(r.applyTo(u.scalarMultiply(v.getNorm())), v.scalarMultiply(u.getNorm()));

        checkAngle(new FieldRotation<Dfp>(u, u.negate()).getAngle(), FastMath.PI);

        try {
            new FieldRotation<Dfp>(u, createVector(0, 0, 0));
            Assert.fail("an exception should have been thrown");
        } catch (MathArithmeticException e) {
            // expected behavior
        }

    }

[620, 'src/test/java', 'org.apache.commons.math3.stat.ranking', 'NaturalRankingTest', 'testNaNsMaximalTiesMinimum', 77, 98]

[620, 'src/test/java', 'org.apache.commons.math3.stat.ranking', 'NaturalRankingTest', 'testNaNsMinimalTiesMaximum', 124, 146]

    @Test
    public void testNaNsMaximalTiesMinimum() {
        NaturalRanking ranking = new NaturalRanking(NaNStrategy.MAXIMAL, TiesStrategy.MINIMUM);
        double[] ranks = ranking.rank(exampleData);
        double[] correctRanks = { 5, 2, 6, 7, 2, 8, 9, 1, 2 };
        TestUtils.assertEquals(correctRanks, ranks, 0d);
        ranks = ranking.rank(tiesFirst);
        correctRanks = new double[] { 1, 1, 4, 3, 5 };
        TestUtils.assertEquals(correctRanks, ranks, 0d);
        ranks = ranking.rank(tiesLast);
        correctRanks = new double[] { 3, 3, 2, 1 };
        TestUtils.assertEquals(correctRanks, ranks, 0d);
        ranks = ranking.rank(multipleNaNs);
        correctRanks = new double[] { 1, 2, 3, 3 };
        TestUtils.assertEquals(correctRanks, ranks, 0d);
        ranks = ranking.rank(multipleTies);
        correctRanks = new double[] { 3, 2, 4, 4, 6, 6, 1 };
        TestUtils.assertEquals(correctRanks, ranks, 0d);
        ranks = ranking.rank(allSame);
        correctRanks = new double[] { 1, 1, 1, 1 };
        TestUtils.assertEquals(correctRanks, ranks, 0d);
    }

    @Test
    public void testNaNsMinimalTiesMaximum() {
        NaturalRanking ranking = new NaturalRanking(NaNStrategy.MINIMAL,
                TiesStrategy.MAXIMUM);
        double[] ranks = ranking.rank(exampleData);
        double[] correctRanks = { 6, 5, 7, 8, 5, 9, 2, 2, 5 };
        TestUtils.assertEquals(correctRanks, ranks, 0d);
        ranks = ranking.rank(tiesFirst);
        correctRanks = new double[] { 2, 2, 4, 3, 5 };
        TestUtils.assertEquals(correctRanks, ranks, 0d);
        ranks = ranking.rank(tiesLast);
        correctRanks = new double[] { 4, 4, 2, 1 };
        TestUtils.assertEquals(correctRanks, ranks, 0d);
        ranks = ranking.rank(multipleNaNs);
        correctRanks = new double[] { 3, 4, 2, 2 };
        TestUtils.assertEquals(correctRanks, ranks, 0d);
        ranks = ranking.rank(multipleTies);
        correctRanks = new double[] { 3, 2, 5, 5, 7, 7, 1 };
        TestUtils.assertEquals(correctRanks, ranks, 0d);
        ranks = ranking.rank(allSame);
        correctRanks = new double[] { 4, 4, 4, 4 };
        TestUtils.assertEquals(correctRanks, ranks, 0d);
    }

[629, 'src/test/java', 'org.apache.commons.math3.fraction', 'BigFractionFormatTest', 'testFormatNegative', 59, 69]

[629, 'src/test/java', 'org.apache.commons.math3.fraction', 'FractionFormatTest', 'testFormatNegative', 57, 67]

    @Test
    public void testFormatNegative() {
        BigFraction c = new BigFraction(-1, 2);
        String expected = "-1 / 2";

        String actual = properFormat.format(c);
        Assert.assertEquals(expected, actual);

        actual = improperFormat.format(c);
        Assert.assertEquals(expected, actual);
    }

    @Test
    public void testFormatNegative() {
        Fraction c = new Fraction(-1, 2);
        String expected = "-1 / 2";

        String actual = properFormat.format(c);
        Assert.assertEquals(expected, actual);

        actual = improperFormat.format(c);
        Assert.assertEquals(expected, actual);
    }

[640, 'src/test/java', 'org.apache.commons.math3.stat.inference', 'ChiSquareTestTest', 'testChiSquareIndependence', 104, 161]

[640, 'src/test/java', 'org.apache.commons.math3.stat.inference', 'TestUtilsTest', 'testChiSquareIndependence', 108, 165]

    @Test
    public void testChiSquareIndependence() {

        // Target values computed using R version 1.8.1

        long[][] counts = { {40, 22, 43}, {91, 21, 28}, {60, 10, 22}};
        Assert.assertEquals( "chi-square test statistic", 22.709027688, testStatistic.chiSquare(counts), 1E-9);
        Assert.assertEquals("chi-square p-value", 0.000144751460134, testStatistic.chiSquareTest(counts), 1E-9);
        Assert.assertTrue("chi-square test reject", testStatistic.chiSquareTest(counts, 0.0002));
        Assert.assertTrue("chi-square test accept", !testStatistic.chiSquareTest(counts, 0.0001));

        long[][] counts2 = {{10, 15}, {30, 40}, {60, 90} };
        Assert.assertEquals( "chi-square test statistic", 0.168965517241, testStatistic.chiSquare(counts2), 1E-9);
        Assert.assertEquals("chi-square p-value",0.918987499852, testStatistic.chiSquareTest(counts2), 1E-9);
        Assert.assertTrue("chi-square test accept", !testStatistic.chiSquareTest(counts2, 0.1));

        // ragged input array
        long[][] counts3 = { {40, 22, 43}, {91, 21, 28}, {60, 10}};
        try {
            testStatistic.chiSquare(counts3);
            Assert.fail("Expecting DimensionMismatchException");
        } catch (DimensionMismatchException ex) {
            // expected
        }

        // insufficient data
        long[][] counts4 = {{40, 22, 43}};
        try {
            testStatistic.chiSquare(counts4);
            Assert.fail("Expecting DimensionMismatchException");
        } catch (DimensionMismatchException ex) {
            // expected
        }
        long[][] counts5 = {{40}, {40}, {30}, {10}};
        try {
            testStatistic.chiSquare(counts5);
            Assert.fail("Expecting DimensionMismatchException");
        } catch (DimensionMismatchException ex) {
            // expected
        }

        // negative counts
        long[][] counts6 = {{10, -2}, {30, 40}, {60, 90} };
        try {
            testStatistic.chiSquare(counts6);
            Assert.fail("Expecting NotPositiveException");
        } catch (NotPositiveException ex) {
            // expected
        }

        // bad alpha
        try {
            testStatistic.chiSquareTest(counts, 0);
            Assert.fail("Expecting OutOfRangeException");
        } catch (OutOfRangeException ex) {
            // expected
        }
    }

    @Test
    public void testChiSquareIndependence() {

        // Target values computed using R version 1.8.1

        long[][] counts = { {40, 22, 43}, {91, 21, 28}, {60, 10, 22}};
        Assert.assertEquals( "chi-square test statistic", 22.709027688, TestUtils.chiSquare(counts), 1E-9);
        Assert.assertEquals("chi-square p-value", 0.000144751460134, TestUtils.chiSquareTest(counts), 1E-9);
        Assert.assertTrue("chi-square test reject", TestUtils.chiSquareTest(counts, 0.0002));
        Assert.assertTrue("chi-square test accept", !TestUtils.chiSquareTest(counts, 0.0001));

        long[][] counts2 = {{10, 15}, {30, 40}, {60, 90} };
        Assert.assertEquals( "chi-square test statistic", 0.168965517241, TestUtils.chiSquare(counts2), 1E-9);
        Assert.assertEquals("chi-square p-value",0.918987499852, TestUtils.chiSquareTest(counts2), 1E-9);
        Assert.assertTrue("chi-square test accept", !TestUtils.chiSquareTest(counts2, 0.1));

        // ragged input array
        long[][] counts3 = { {40, 22, 43}, {91, 21, 28}, {60, 10}};
        try {
            TestUtils.chiSquare(counts3);
            Assert.fail("Expecting DimensionMismatchException");
        } catch (DimensionMismatchException ex) {
            // expected
        }

        // insufficient data
        long[][] counts4 = {{40, 22, 43}};
        try {
            TestUtils.chiSquare(counts4);
            Assert.fail("Expecting DimensionMismatchException");
        } catch (DimensionMismatchException ex) {
            // expected
        }
        long[][] counts5 = {{40}, {40}, {30}, {10}};
        try {
            TestUtils.chiSquare(counts5);
            Assert.fail("Expecting DimensionMismatchException");
        } catch (DimensionMismatchException ex) {
            // expected
        }

        // negative counts
        long[][] counts6 = {{10, -2}, {30, 40}, {60, 90} };
        try {
            TestUtils.chiSquare(counts6);
            Assert.fail("Expecting NotPositiveException");
        } catch (NotPositiveException ex) {
            // expected
        }

        // bad alpha
        try {
            TestUtils.chiSquareTest(counts, 0);
            Assert.fail("Expecting OutOfRangeException");
        } catch (OutOfRangeException ex) {
            // expected
        }
    }

[649, 'src/test/java', 'org.apache.commons.math3.distribution', 'LogNormalDistributionTest', 'testExtremeValues', 207, 224]

[649, 'src/test/java', 'org.apache.commons.math3.distribution', 'ParetoDistributionTest', 'testExtremeValues', 173, 190]

    @Test
    public void testExtremeValues() {
        LogNormalDistribution d = new LogNormalDistribution(0, 1);
        for (int i = 0; i < 1e5; i++) { // make sure no convergence exception
            double upperTail = d.cumulativeProbability(i);
            if (i <= 72) { // make sure not top-coded
                Assert.assertTrue(upperTail < 1.0d);
            }
            else { // make sure top coding not reversed
                Assert.assertTrue(upperTail > 0.99999);
            }
        }

        Assert.assertEquals(d.cumulativeProbability(Double.MAX_VALUE), 1, 0);
        Assert.assertEquals(d.cumulativeProbability(-Double.MAX_VALUE), 0, 0);
        Assert.assertEquals(d.cumulativeProbability(Double.POSITIVE_INFINITY), 1, 0);
        Assert.assertEquals(d.cumulativeProbability(Double.NEGATIVE_INFINITY), 0, 0);
    }

    @Test
    public void testExtremeValues() {
        ParetoDistribution d = new ParetoDistribution(1, 1);
        for (int i = 0; i < 1e5; i++) { // make sure no convergence exception
            double upperTail = d.cumulativeProbability(i);
            if (i <= 1000) { // make sure not top-coded
                Assert.assertTrue(upperTail < 1.0d);
            }
            else { // make sure top coding not reversed
                Assert.assertTrue(upperTail > 0.999);
            }
        }

        Assert.assertEquals(d.cumulativeProbability(Double.MAX_VALUE), 1, 0);
        Assert.assertEquals(d.cumulativeProbability(-Double.MAX_VALUE), 0, 0);
        Assert.assertEquals(d.cumulativeProbability(Double.POSITIVE_INFINITY), 1, 0);
        Assert.assertEquals(d.cumulativeProbability(Double.NEGATIVE_INFINITY), 0, 0);
    }

[693, 'src/test/java', 'org.apache.commons.math3.analysis.interpolation', 'DividedDifferenceInterpolatorTest', 'testSinFunction', 48, 73]

[693, 'src/test/java', 'org.apache.commons.math3.analysis.interpolation', 'NevilleInterpolatorTest', 'testSinFunction', 48, 73]

    @Test
    public void testSinFunction() {
        UnivariateFunction f = new Sin();
        UnivariateInterpolator interpolator = new DividedDifferenceInterpolator();
        double x[], y[], z, expected, result, tolerance;

        // 6 interpolating points on interval [0, 2*PI]
        int n = 6;
        double min = 0.0, max = 2 * FastMath.PI;
        x = new double[n];
        y = new double[n];
        for (int i = 0; i < n; i++) {
            x[i] = min + i * (max - min) / n;
            y[i] = f.value(x[i]);
        }
        double derivativebound = 1.0;
        UnivariateFunction p = interpolator.interpolate(x, y);

        z = FastMath.PI / 4; expected = f.value(z); result = p.value(z);
        tolerance = FastMath.abs(derivativebound * partialerror(x, z));
        Assert.assertEquals(expected, result, tolerance);

        z = FastMath.PI * 1.5; expected = f.value(z); result = p.value(z);
        tolerance = FastMath.abs(derivativebound * partialerror(x, z));
        Assert.assertEquals(expected, result, tolerance);
    }

    @Test
    public void testSinFunction() {
        UnivariateFunction f = new Sin();
        UnivariateInterpolator interpolator = new NevilleInterpolator();
        double x[], y[], z, expected, result, tolerance;

        // 6 interpolating points on interval [0, 2*PI]
        int n = 6;
        double min = 0.0, max = 2 * FastMath.PI;
        x = new double[n];
        y = new double[n];
        for (int i = 0; i < n; i++) {
            x[i] = min + i * (max - min) / n;
            y[i] = f.value(x[i]);
        }
        double derivativebound = 1.0;
        UnivariateFunction p = interpolator.interpolate(x, y);

        z = FastMath.PI / 4; expected = f.value(z); result = p.value(z);
        tolerance = FastMath.abs(derivativebound * partialerror(x, z));
        Assert.assertEquals(expected, result, tolerance);

        z = FastMath.PI * 1.5; expected = f.value(z); result = p.value(z);
        tolerance = FastMath.abs(derivativebound * partialerror(x, z));
        Assert.assertEquals(expected, result, tolerance);
    }

[721, 'src/test/java', 'org.apache.commons.math3.genetics', 'UniformCrossoverTest', 'testCrossoverInvalidFixedLengthChromosomeFirst', 121, 135]

[721, 'src/test/java', 'org.apache.commons.math3.genetics', 'UniformCrossoverTest', 'testCrossoverInvalidFixedLengthChromosomeSecond', 137, 151]

    @Test(expected = MathIllegalArgumentException.class)
    public void testCrossoverInvalidFixedLengthChromosomeFirst() {
        @SuppressWarnings("boxing")
        final Integer[] p1 = new Integer[] {1,0,1,0,0,1,0,1,1};
        final BinaryChromosome p1c = new DummyBinaryChromosome(p1);
        final Chromosome p2c = new Chromosome() {
            public double fitness() {
                // Not important
                return 0;
            }
        };

        final CrossoverPolicy cp = new UniformCrossover<Integer>(0.5d);
        cp.crossover(p1c, p2c);
    }

    @Test(expected = MathIllegalArgumentException.class)
    public void testCrossoverInvalidFixedLengthChromosomeSecond() {
        @SuppressWarnings("boxing")
        final Integer[] p1 = new Integer[] {1,0,1,0,0,1,0,1,1};
        final BinaryChromosome p2c = new DummyBinaryChromosome(p1);
        final Chromosome p1c = new Chromosome() {
            public double fitness() {
                // Not important
                return 0;
            }
        };

        final CrossoverPolicy cp = new UniformCrossover<Integer>(0.5d);
        cp.crossover(p1c, p2c);
    }

[761, 'src/test/java', 'org.apache.commons.math3.linear', 'MatrixUtilsTest', 'testCreateRowRealMatrix', 119, 135]

[761, 'src/test/java', 'org.apache.commons.math3.linear', 'MatrixUtilsTest', 'testCreateColumnRealMatrix', 157, 173]

    @Test
    public void testCreateRowRealMatrix() {
        Assert.assertEquals(MatrixUtils.createRowRealMatrix(row),
                     new BlockRealMatrix(rowMatrix));
        try {
            MatrixUtils.createRowRealMatrix(new double[] {});  // empty
            Assert.fail("Expecting NotStrictlyPositiveException");
        } catch (NotStrictlyPositiveException ex) {
            // expected
        }
        try {
            MatrixUtils.createRowRealMatrix(null);  // null
            Assert.fail("Expecting NullArgumentException");
        } catch (NullArgumentException ex) {
            // expected
        }
    }

    @Test
    public void testCreateColumnRealMatrix() {
        Assert.assertEquals(MatrixUtils.createColumnRealMatrix(col),
                     new BlockRealMatrix(colMatrix));
        try {
            MatrixUtils.createColumnRealMatrix(new double[] {});  // empty
            Assert.fail("Expecting MathIllegalArgumentException");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }
        try {
            MatrixUtils.createColumnRealMatrix(null);  // null
            Assert.fail("Expecting NullArgumentException");
        } catch (NullArgumentException ex) {
            // expected
        }
    }

[786, 'src/test/java', 'org.apache.commons.math3.complex', 'ComplexTest', 'testAddInf', 118, 128]

[786, 'src/test/java', 'org.apache.commons.math3.complex', 'ComplexTest', 'testSubtractInf', 453, 463]

    @Test
    public void testAddInf() {
        Complex x = new Complex(1, 1);
        Complex z = new Complex(inf, 0);
        Complex w = x.add(z);
        Assert.assertEquals(w.getImaginary(), 1, 0);
        Assert.assertEquals(inf, w.getReal(), 0);

        x = new Complex(neginf, 0);
        Assert.assertTrue(Double.isNaN(x.add(z).getReal()));
    }

    @Test
    public void testSubtractInf() {
        Complex x = new Complex(1, 1);
        Complex z = new Complex(neginf, 0);
        Complex w = x.subtract(z);
        Assert.assertEquals(w.getImaginary(), 1, 0);
        Assert.assertEquals(inf, w.getReal(), 0);

        x = new Complex(neginf, 0);
        Assert.assertTrue(Double.isNaN(x.subtract(z).getReal()));
    }

[788, 'src/test/java', 'org.apache.commons.math3.geometry.euclidean.threed', 'FieldRotationDSTest', 'testIssue639', 1187, 1200]

[788, 'src/test/java', 'org.apache.commons.math3.geometry.euclidean.threed', 'FieldRotationDfpTest', 'testIssue639', 952, 965]

    @Test
    public void testIssue639() throws MathArithmeticException{
        FieldVector3D<DerivativeStructure> u1 = createVector(-1321008684645961.0 /  268435456.0,
                                   -5774608829631843.0 /  268435456.0,
                                   -3822921525525679.0 / 4294967296.0);
        FieldVector3D<DerivativeStructure> u2 =createVector( -5712344449280879.0 /    2097152.0,
                                   -2275058564560979.0 /    1048576.0,
                                   4423475992255071.0 /      65536.0);
        FieldRotation<DerivativeStructure> rot = new FieldRotation<DerivativeStructure>(u1, u2, createVector(1, 0, 0),createVector(0, 0, 1));
        Assert.assertEquals( 0.6228370359608200639829222, rot.getQ0().getReal(), 1.0e-15);
        Assert.assertEquals( 0.0257707621456498790029987, rot.getQ1().getReal(), 1.0e-15);
        Assert.assertEquals(-0.0000000002503012255839931, rot.getQ2().getReal(), 1.0e-15);
        Assert.assertEquals(-0.7819270390861109450724902, rot.getQ3().getReal(), 1.0e-15);
    }

    @Test
    public void testIssue639() throws MathArithmeticException{
        FieldVector3D<Dfp> u1 = createVector(-1321008684645961.0 /  268435456.0,
                                   -5774608829631843.0 /  268435456.0,
                                   -3822921525525679.0 / 4294967296.0);
        FieldVector3D<Dfp> u2 =createVector( -5712344449280879.0 /    2097152.0,
                                   -2275058564560979.0 /    1048576.0,
                                   4423475992255071.0 /      65536.0);
        FieldRotation<Dfp> rot = new FieldRotation<Dfp>(u1, u2, createVector(1, 0, 0),createVector(0, 0, 1));
        Assert.assertEquals( 0.6228370359608200639829222, rot.getQ0().getReal(), 1.0e-15);
        Assert.assertEquals( 0.0257707621456498790029987, rot.getQ1().getReal(), 1.0e-15);
        Assert.assertEquals(-0.0000000002503012255839931, rot.getQ2().getReal(), 1.0e-15);
        Assert.assertEquals(-0.7819270390861109450724902, rot.getQ3().getReal(), 1.0e-15);
    }

[806, 'src/test/java', 'org.apache.commons.math3.linear', 'ConjugateGradientTest', 'testUnpreconditionedInPlaceSolutionWithInitialGuess', 101, 125]

[806, 'src/test/java', 'org.apache.commons.math3.linear', 'SymmLQTest', 'testUnpreconditionedInPlaceSolutionWithInitialGuess', 267, 291]

    @Test
    public void testUnpreconditionedInPlaceSolutionWithInitialGuess() {
        final int n = 5;
        final int maxIterations = 100;
        final RealLinearOperator a = new HilbertMatrix(n);
        final InverseHilbertMatrix ainv = new InverseHilbertMatrix(n);
        final IterativeLinearSolver solver;
        solver = new ConjugateGradient(maxIterations, 1E-10, true);
        final RealVector b = new ArrayRealVector(n);
        for (int j = 0; j < n; j++) {
            b.set(0.);
            b.setEntry(j, 1.);
            final RealVector x0 = new ArrayRealVector(n);
            x0.set(1.);
            final RealVector x = solver.solveInPlace(a, b, x0);
            Assert.assertSame("x should be a reference to x0", x0, x);
            for (int i = 0; i < n; i++) {
                final double actual = x.getEntry(i);
                final double expected = ainv.getEntry(i, j);
                final double delta = 1E-10 * FastMath.abs(expected);
                final String msg = String.format("entry[%d][%d)", i, j);
                Assert.assertEquals(msg, expected, actual, delta);
            }
        }
    }

    @Test
    public void testUnpreconditionedInPlaceSolutionWithInitialGuess() {
        final int n = 5;
        final int maxIterations = 100;
        final RealLinearOperator a = new HilbertMatrix(n);
        final InverseHilbertMatrix ainv = new InverseHilbertMatrix(n);
        final IterativeLinearSolver solver;
        solver = new SymmLQ(maxIterations, 1E-10, true);
        final RealVector b = new ArrayRealVector(n);
        for (int j = 0; j < n; j++) {
            b.set(0.);
            b.setEntry(j, 1.);
            final RealVector x0 = new ArrayRealVector(n);
            x0.set(1.);
            final RealVector x = solver.solveInPlace(a, b, x0);
            Assert.assertSame("x should be a reference to x0", x0, x);
            for (int i = 0; i < n; i++) {
                final double actual = x.getEntry(i);
                final double expected = ainv.getEntry(i, j);
                final double delta = 1E-6 * FastMath.abs(expected);
                final String msg = String.format("entry[%d][%d)", i, j);
                Assert.assertEquals(msg, expected, actual, delta);
            }
        }
    }

[809, 'src/test/java', 'org.apache.commons.math3.linear', 'BlockFieldMatrixTest', 'testGetSetRowVectorLarge', 911, 929]

[809, 'src/test/java', 'org.apache.commons.math3.linear', 'BlockFieldMatrixTest', 'testGetSetColumnVectorLarge', 973, 991]

    @Test
    public void testGetSetRowVectorLarge() {
        int n = 3 * BlockFieldMatrix.BLOCK_SIZE;
        FieldMatrix<Fraction> m = new BlockFieldMatrix<Fraction>(FractionField.getInstance(), n, n);
        FieldVector<Fraction> sub = new ArrayFieldVector<Fraction>(n, new Fraction(1));

        m.setRowVector(2, sub);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                if (i != 2) {
                    Assert.assertEquals(new Fraction(0), m.getEntry(i, j));
                } else {
                    Assert.assertEquals(new Fraction(1), m.getEntry(i, j));
                }
            }
        }
        Assert.assertEquals(sub, m.getRowVector(2));

    }

    @Test
    public void testGetSetColumnVectorLarge() {
        int n = 3 * BlockFieldMatrix.BLOCK_SIZE;
        FieldMatrix<Fraction> m = new BlockFieldMatrix<Fraction>(FractionField.getInstance(), n, n);
        FieldVector<Fraction> sub = new ArrayFieldVector<Fraction>(n, new Fraction(1));

        m.setColumnVector(2, sub);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                if (j != 2) {
                    Assert.assertEquals(new Fraction(0), m.getEntry(i, j));
                } else {
                    Assert.assertEquals(new Fraction(1), m.getEntry(i, j));
                }
            }
        }
        Assert.assertEquals(sub, m.getColumnVector(2));

    }

[831, 'src/test/java', 'org.apache.commons.math3.linear', 'BlockFieldMatrixTest', 'testGetSetRowMatrixLarge', 783, 803]

[831, 'src/test/java', 'org.apache.commons.math3.linear', 'BlockFieldMatrixTest', 'testGetSetColumnMatrixLarge', 847, 867]

    @Test
    public void testGetSetRowMatrixLarge() {
        int n = 3 * BlockFieldMatrix.BLOCK_SIZE;
        FieldMatrix<Fraction> m =
            new BlockFieldMatrix<Fraction>(FractionField.getInstance(), n, n);
        FieldMatrix<Fraction> sub =
            new BlockFieldMatrix<Fraction>(FractionField.getInstance(), 1, n).scalarAdd(new Fraction(1));

        m.setRowMatrix(2, sub);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                if (i != 2) {
                    Assert.assertEquals(new Fraction(0), m.getEntry(i, j));
                } else {
                    Assert.assertEquals(new Fraction(1), m.getEntry(i, j));
                }
            }
        }
        Assert.assertEquals(sub, m.getRowMatrix(2));

    }

    @Test
    public void testGetSetColumnMatrixLarge() {
        int n = 3 * BlockFieldMatrix.BLOCK_SIZE;
        FieldMatrix<Fraction> m =
            new BlockFieldMatrix<Fraction>(FractionField.getInstance(), n, n);
        FieldMatrix<Fraction> sub =
            new BlockFieldMatrix<Fraction>(FractionField.getInstance(), n, 1).scalarAdd(new Fraction(1));

        m.setColumnMatrix(2, sub);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                if (j != 2) {
                    Assert.assertEquals(new Fraction(0), m.getEntry(i, j));
                } else {
                    Assert.assertEquals(new Fraction(1), m.getEntry(i, j));
                }
            }
        }
        Assert.assertEquals(sub, m.getColumnMatrix(2));

    }

[895, 'src/test/java', 'org.apache.commons.math3.ode.nonstiff', 'AdamsBashforthIntegratorTest', 'exceedMaxEvaluations', 111, 125]

[895, 'src/test/java', 'org.apache.commons.math3.ode.nonstiff', 'AdamsMoultonIntegratorTest', 'exceedMaxEvaluations', 116, 132]

    @Test(expected = MaxCountExceededException.class)
    public void exceedMaxEvaluations() throws DimensionMismatchException, NumberIsTooSmallException, MaxCountExceededException, NoBracketingException {

        TestProblem1 pb  = new TestProblem1();
        double range = pb.getFinalTime() - pb.getInitialTime();

        AdamsBashforthIntegrator integ = new AdamsBashforthIntegrator(2, 0, range, 1.0e-12, 1.0e-12);
        TestProblemHandler handler = new TestProblemHandler(pb, integ);
        integ.addStepHandler(handler);
        integ.setMaxEvaluations(650);
        integ.integrate(pb,
                        pb.getInitialTime(), pb.getInitialState(),
                        pb.getFinalTime(), new double[pb.getDimension()]);

    }

    @Test(expected = MaxCountExceededException.class)
    public void exceedMaxEvaluations()
            throws DimensionMismatchException, NumberIsTooSmallException,
            MaxCountExceededException, NoBracketingException {

        TestProblem1 pb  = new TestProblem1();
        double range = pb.getFinalTime() - pb.getInitialTime();

        AdamsMoultonIntegrator integ = new AdamsMoultonIntegrator(2, 0, range, 1.0e-12, 1.0e-12);
        TestProblemHandler handler = new TestProblemHandler(pb, integ);
        integ.addStepHandler(handler);
        integ.setMaxEvaluations(650);
        integ.integrate(pb,
                        pb.getInitialTime(), pb.getInitialState(),
                        pb.getFinalTime(), new double[pb.getDimension()]);

    }

[896, 'src/test/java', 'org.apache.commons.math3.analysis.differentiation', 'SparseGradientTest', 'testSinhDefinition', 701, 710]

[896, 'src/test/java', 'org.apache.commons.math3.analysis.differentiation', 'SparseGradientTest', 'testCoshDefinition', 712, 721]

    @Test
    public void testSinhDefinition() {
        for (double x = 0.1; x < 1.2; x += 0.001) {
            SparseGradient sgX = SparseGradient.createVariable(0, x);
            SparseGradient sinh1 = sgX.exp().subtract(sgX.exp().reciprocal()).multiply(0.5);
            SparseGradient sinh2 = sgX.sinh();
            SparseGradient zero = sinh1.subtract(sinh2);
            checkF0F1(zero, 0.0, 0.0);
        }
    }

    @Test
    public void testCoshDefinition() {
        for (double x = 0.1; x < 1.2; x += 0.001) {
            SparseGradient sgX = SparseGradient.createVariable(0, x);
            SparseGradient cosh1 = sgX.exp().add(sgX.exp().reciprocal()).multiply(0.5);
            SparseGradient cosh2 = sgX.cosh();
            SparseGradient zero = cosh1.subtract(cosh2);
            checkF0F1(zero, 0.0, 0.0);
        }
    }

[905, 'src/test/java', 'org.apache.commons.math3.linear', 'BlockFieldMatrixTest', 'testSerial', 1329, 1333]

[905, 'src/test/java', 'org.apache.commons.math3.linear', 'FieldMatrixImplTest', 'testSerial', 1056, 1060]

    @Test
    public void testSerial()  {
        BlockFieldMatrix<Fraction> m = new BlockFieldMatrix<Fraction>(testData);
        Assert.assertEquals(m,TestUtils.serializeAndRecover(m));
    }

    @Test
    public void testSerial()  {
        Array2DRowFieldMatrix<Fraction> m = new Array2DRowFieldMatrix<Fraction>(testData);
        Assert.assertEquals(m,TestUtils.serializeAndRecover(m));
    }

[908, 'src/test/java', 'org.apache.commons.math3.complex', 'ComplexTest', 'testSin', 996, 1001]

[908, 'src/test/java', 'org.apache.commons.math3.complex', 'ComplexTest', 'testSqrt1z', 1111, 1116]

    @Test
    public void testSin() {
        Complex z = new Complex(3, 4);
        Complex expected = new Complex(3.853738, -27.01681);
        TestUtils.assertEquals(expected, z.sin(), 1.0e-5);
    }

    @Test
    public void testSqrt1z() {
        Complex z = new Complex(3, 4);
        Complex expected = new Complex(4.08033, -2.94094);
        TestUtils.assertEquals(expected, z.sqrt1z(), 1.0e-5);
    }

[914, 'src/test/java', 'org.apache.commons.math3.random', 'HaltonSequenceGeneratorTest', 'test3DReference', 62, 69]

[914, 'src/test/java', 'org.apache.commons.math3.random', 'SobolSequenceGeneratorTest', 'test3DReference', 49, 56]

    @Test
    public void test3DReference() {
        for (int i = 0; i < referenceValues.length; i++) {
            double[] result = generator.nextVector();
            Assert.assertArrayEquals(referenceValues[i], result, 1e-3);
            Assert.assertEquals(i + 1, generator.getNextIndex());
        }
    }

    @Test
    public void test3DReference() {
        for (int i = 0; i < referenceValues.length; i++) {
            double[] result = generator.nextVector();
            Assert.assertArrayEquals(referenceValues[i], result, 1e-6);
            Assert.assertEquals(i + 1, generator.getNextIndex());
        }
    }

[938, 'src/test/java', 'org.apache.commons.math3.linear', 'Array2DRowRealMatrixTest', 'testDimensions', 101, 111]

[938, 'src/test/java', 'org.apache.commons.math3.linear', 'BlockRealMatrixTest', 'testDimensions', 101, 111]

    @Test
    public void testDimensions() {
        Array2DRowRealMatrix m = new Array2DRowRealMatrix(testData);
        Array2DRowRealMatrix m2 = new Array2DRowRealMatrix(testData2);
        Assert.assertEquals("testData row dimension",3,m.getRowDimension());
        Assert.assertEquals("testData column dimension",3,m.getColumnDimension());
        Assert.assertTrue("testData is square",m.isSquare());
        Assert.assertEquals("testData2 row dimension",m2.getRowDimension(),2);
        Assert.assertEquals("testData2 column dimension",m2.getColumnDimension(),3);
        Assert.assertTrue("testData2 is not square",!m2.isSquare());
    }

    @Test
    public void testDimensions() {
        BlockRealMatrix m = new BlockRealMatrix(testData);
        BlockRealMatrix m2 = new BlockRealMatrix(testData2);
        Assert.assertEquals("testData row dimension",3,m.getRowDimension());
        Assert.assertEquals("testData column dimension",3,m.getColumnDimension());
        Assert.assertTrue("testData is square",m.isSquare());
        Assert.assertEquals("testData2 row dimension",m2.getRowDimension(),2);
        Assert.assertEquals("testData2 column dimension",m2.getColumnDimension(),3);
        Assert.assertTrue("testData2 is not square",!m2.isSquare());
    }

[939, 'src/test/java', 'org.apache.commons.math3.stat.descriptive.rank', 'PSquarePercentileTest', 'testMarkersWithLowerIndex', 272, 283]

[939, 'src/test/java', 'org.apache.commons.math3.stat.descriptive.rank', 'PSquarePercentileTest', 'testMarkersWithHigherIndex', 285, 296]

    @Test(expected = OutOfRangeException.class)
    public void testMarkersWithLowerIndex() {
        PSquareMarkers mThat =
                PSquarePercentile.newMarkers(
                        Arrays.asList(new Double[] { 95.1772, 95.1567, 95.1937,
                                95.1959, 95.1442, 95.0610, 95.1591, 95.1195,
                                95.1772, 95.0925, 95.1990, 95.1682 }), 0.50);
        for (int i = 0; i < testArray.length; i++) {
            mThat.processDataPoint(testArray[i]);
        }
        mThat.estimate(0);
    }

    @Test(expected = OutOfRangeException.class)
    public void testMarkersWithHigherIndex() {
        PSquareMarkers mThat =
                PSquarePercentile.newMarkers(
                        Arrays.asList(new Double[] { 95.1772, 95.1567, 95.1937,
                                95.1959, 95.1442, 95.0610, 95.1591, 95.1195,
                                95.1772, 95.0925, 95.1990, 95.1682 }), 0.50);
        for (int i = 0; i < testArray.length; i++) {
            mThat.processDataPoint(testArray[i]);
        }
        mThat.estimate(6);
    }

[962, 'src/test/java', 'org.apache.commons.math3.analysis.function', 'GaussianTest', 'testParametricUsage3', 107, 111]

[962, 'src/test/java', 'org.apache.commons.math3.analysis.function', 'GaussianTest', 'testParametricUsage6', 125, 129]

    @Test(expected=NotStrictlyPositiveException.class)
    public void testParametricUsage3() {
        final Gaussian.Parametric g = new Gaussian.Parametric();
        g.value(0, new double[] {0, 1, 0});
    }

    @Test(expected=NotStrictlyPositiveException.class)
    public void testParametricUsage6() {
        final Gaussian.Parametric g = new Gaussian.Parametric();
        g.gradient(0, new double[] {0, 1, 0});
    }

[1005, 'src/test/java', 'org.apache.commons.math3.optim.linear', 'SimplexSolverTest', 'testMath272', 219, 235]

[1005, 'src/test/java', 'org.apache.commons.math3.optim.linear', 'SimplexSolverTest', 'testEpsilon', 592, 608]

    @Test
    public void testMath272() {
        LinearObjectiveFunction f = new LinearObjectiveFunction(new double[] { 2, 2, 1 }, 0);
        Collection<LinearConstraint> constraints = new ArrayList<LinearConstraint>();
        constraints.add(new LinearConstraint(new double[] { 1, 1, 0 }, Relationship.GEQ,  1));
        constraints.add(new LinearConstraint(new double[] { 1, 0, 1 }, Relationship.GEQ,  1));
        constraints.add(new LinearConstraint(new double[] { 0, 1, 0 }, Relationship.GEQ,  1));

        SimplexSolver solver = new SimplexSolver();
        PointValuePair solution = solver.optimize(DEFAULT_MAX_ITER, f, new LinearConstraintSet(constraints),
                                                  GoalType.MINIMIZE, new NonNegativeConstraint(true));

        Assert.assertEquals(0.0, solution.getPoint()[0], .0000001);
        Assert.assertEquals(1.0, solution.getPoint()[1], .0000001);
        Assert.assertEquals(1.0, solution.getPoint()[2], .0000001);
        Assert.assertEquals(3.0, solution.getValue(), .0000001);
    }

    @Test
    public void testEpsilon() {
      LinearObjectiveFunction f =
          new LinearObjectiveFunction(new double[] { 10, 5, 1 }, 0);
      Collection<LinearConstraint> constraints = new ArrayList<LinearConstraint>();
      constraints.add(new LinearConstraint(new double[] {  9, 8, 0 }, Relationship.EQ,  17));
      constraints.add(new LinearConstraint(new double[] {  0, 7, 8 }, Relationship.LEQ,  7));
      constraints.add(new LinearConstraint(new double[] { 10, 0, 2 }, Relationship.LEQ, 10));

      SimplexSolver solver = new SimplexSolver();
      PointValuePair solution = solver.optimize(DEFAULT_MAX_ITER, f, new LinearConstraintSet(constraints),
                                                GoalType.MAXIMIZE, new NonNegativeConstraint(false));
      Assert.assertEquals(1.0, solution.getPoint()[0], 0.0);
      Assert.assertEquals(1.0, solution.getPoint()[1], 0.0);
      Assert.assertEquals(0.0, solution.getPoint()[2], 0.0);
      Assert.assertEquals(15.0, solution.getValue(), 0.0);
  }

[1018, 'src/test/java', 'org.apache.commons.math3.ode.nonstiff', 'ClassicalRungeKuttaIntegratorTest', 'testMissedEndEvent', 44, 101]

[1018, 'src/test/java', 'org.apache.commons.math3.ode.nonstiff', 'LutherIntegratorTest', 'testMissedEndEvent', 44, 101]

  @Test
  public void testMissedEndEvent()
      throws DimensionMismatchException, NumberIsTooSmallException,
             MaxCountExceededException, NoBracketingException {
      final double   t0     = 1878250320.0000029;
      final double   tEvent = 1878250379.9999986;
      final double[] k      = { 1.0e-4, 1.0e-5, 1.0e-6 };
      FirstOrderDifferentialEquations ode = new FirstOrderDifferentialEquations() {

          public int getDimension() {
              return k.length;
          }

          public void computeDerivatives(double t, double[] y, double[] yDot) {
              for (int i = 0; i < y.length; ++i) {
                  yDot[i] = k[i] * y[i];
              }
          }
      };

      ClassicalRungeKuttaIntegrator integrator = new ClassicalRungeKuttaIntegrator(60.0);

      double[] y0   = new double[k.length];
      for (int i = 0; i < y0.length; ++i) {
          y0[i] = i + 1;
      }
      double[] y    = new double[k.length];

      double finalT = integrator.integrate(ode, t0, y0, tEvent, y);
      Assert.assertEquals(tEvent, finalT, 5.0e-6);
      for (int i = 0; i < y.length; ++i) {
          Assert.assertEquals(y0[i] * FastMath.exp(k[i] * (finalT - t0)), y[i], 1.0e-9);
      }

      integrator.addEventHandler(new EventHandler() {

          public void init(double t0, double[] y0, double t) {
          }

          public void resetState(double t, double[] y) {
          }

          public double g(double t, double[] y) {
              return t - tEvent;
          }

          public Action eventOccurred(double t, double[] y, boolean increasing) {
              Assert.assertEquals(tEvent, t, 5.0e-6);
              return Action.CONTINUE;
          }
      }, Double.POSITIVE_INFINITY, 1.0e-20, 100);
      finalT = integrator.integrate(ode, t0, y0, tEvent + 120, y);
      Assert.assertEquals(tEvent + 120, finalT, 5.0e-6);
      for (int i = 0; i < y.length; ++i) {
          Assert.assertEquals(y0[i] * FastMath.exp(k[i] * (finalT - t0)), y[i], 1.0e-9);
      }

  }

    @Test
    public void testMissedEndEvent()
            throws DimensionMismatchException, NumberIsTooSmallException,
            MaxCountExceededException, NoBracketingException {
        final double   t0     = 1878250320.0000029;
        final double   tEvent = 1878250379.9999986;
        final double[] k      = { 1.0e-4, 1.0e-5, 1.0e-6 };
        FirstOrderDifferentialEquations ode = new FirstOrderDifferentialEquations() {

            public int getDimension() {
                return k.length;
            }

            public void computeDerivatives(double t, double[] y, double[] yDot) {
                for (int i = 0; i < y.length; ++i) {
                    yDot[i] = k[i] * y[i];
                }
            }
        };

        LutherIntegrator integrator = new LutherIntegrator(60.0);

        double[] y0   = new double[k.length];
        for (int i = 0; i < y0.length; ++i) {
            y0[i] = i + 1;
        }
        double[] y    = new double[k.length];

        double finalT = integrator.integrate(ode, t0, y0, tEvent, y);
        Assert.assertEquals(tEvent, finalT, 1.0e-15);
        for (int i = 0; i < y.length; ++i) {
            Assert.assertEquals(y0[i] * FastMath.exp(k[i] * (finalT - t0)), y[i], 1.0e-15);
        }

        integrator.addEventHandler(new EventHandler() {

            public void init(double t0, double[] y0, double t) {
            }

            public void resetState(double t, double[] y) {
            }

            public double g(double t, double[] y) {
                return t - tEvent;
            }

            public Action eventOccurred(double t, double[] y, boolean increasing) {
                Assert.assertEquals(tEvent, t, 1.0e-15);
                return Action.CONTINUE;
            }
        }, Double.POSITIVE_INFINITY, 1.0e-20, 100);
        finalT = integrator.integrate(ode, t0, y0, tEvent + 120, y);
        Assert.assertEquals(tEvent + 120, finalT, 1.0e-15);
        for (int i = 0; i < y.length; ++i) {
            Assert.assertEquals(y0[i] * FastMath.exp(k[i] * (finalT - t0)), y[i], 1.0e-15);
        }

    }

[1031, 'src/test/java', 'org.apache.commons.math3.linear', 'Array2DRowRealMatrixTest', 'testSetColumn', 854, 873]

[1031, 'src/test/java', 'org.apache.commons.math3.linear', 'BlockRealMatrixTest', 'testSetColumn', 993, 1012]

    @Test
    public void testSetColumn() {
        RealMatrix m = new Array2DRowRealMatrix(subTestData);
        double[] mColumn3 = columnToArray(subColumn3);
        Assert.assertTrue(mColumn3[0] != m.getColumn(1)[0]);
        m.setColumn(1, mColumn3);
        checkArrays(mColumn3, m.getColumn(1));
        try {
            m.setColumn(-1, mColumn3);
            Assert.fail("Expecting OutOfRangeException");
        } catch (OutOfRangeException ex) {
            // expected
        }
        try {
            m.setColumn(0, new double[5]);
            Assert.fail("Expecting MatrixDimensionMismatchException");
        } catch (MatrixDimensionMismatchException ex) {
            // expected
        }
    }

    @Test
    public void testSetColumn() {
        RealMatrix m = new BlockRealMatrix(subTestData);
        double[] mColumn3 = columnToArray(subColumn3);
        Assert.assertTrue(mColumn3[0] != m.getColumn(1)[0]);
        m.setColumn(1, mColumn3);
        checkArrays(mColumn3, m.getColumn(1));
        try {
            m.setColumn(-1, mColumn3);
            Assert.fail("Expecting OutOfRangeException");
        } catch (OutOfRangeException ex) {
            // expected
        }
        try {
            m.setColumn(0, new double[5]);
            Assert.fail("Expecting MatrixDimensionMismatchException");
        } catch (MatrixDimensionMismatchException ex) {
            // expected
        }
    }

[1037, 'src/test/java', 'org.apache.commons.math3.distribution', 'TriangularDistributionTest', 'testGetUpperBound', 146, 150]

[1037, 'src/test/java', 'org.apache.commons.math3.distribution', 'UniformRealDistributionTest', 'testGetUpperBound', 78, 82]

    @Test
    public void testGetUpperBound() {
        TriangularDistribution distribution = makeDistribution();
        Assert.assertEquals(12.0, distribution.getSupportUpperBound(), 0);
    }

    @Test
    public void testGetUpperBound() {
        UniformRealDistribution distribution = makeDistribution();
        Assert.assertEquals(1.25, distribution.getSupportUpperBound(), 0);
    }

[1043, 'src/test/java', 'org.apache.commons.math3.optim.nonlinear.vector.jacobian', 'GaussNewtonOptimizerTest', 'testConstraintsUnsupported', 104, 112]

[1043, 'src/test/java', 'org.apache.commons.math3.optim.nonlinear.vector.jacobian', 'LevenbergMarquardtOptimizerTest', 'testConstraintsUnsupported', 113, 121]

    @Test(expected=MathUnsupportedOperationException.class)
    public void testConstraintsUnsupported() {
        createOptimizer().optimize(new MaxEval(100),
                                   new Target(new double[] { 2 }),
                                   new Weight(new double[] { 1 }),
                                   new InitialGuess(new double[] { 1, 2 }),
                                   new SimpleBounds(new double[] { -10, 0 },
                                                    new double[] { 20, 30 }));
    }

    @Test(expected=MathUnsupportedOperationException.class)
    public void testConstraintsUnsupported() {
        createOptimizer().optimize(new MaxEval(100),
                                   new Target(new double[] { 2 }),
                                   new Weight(new double[] { 1 }),
                                   new InitialGuess(new double[] { 1, 2 }),
                                   new SimpleBounds(new double[] { -10, 0 },
                                                    new double[] { 20, 30 }));
    }

[1055, 'src/test/java', 'org.apache.commons.math3.geometry.euclidean.threed', 'RotationTest', 'testRevertVectorOperator', 162, 173]

[1055, 'src/test/java', 'org.apache.commons.math3.geometry.euclidean.threed', 'RotationTest', 'testRevertFrameTransform', 175, 186]

  @Test
  public void testRevertVectorOperator() {
    Rotation r = new Rotation(0.001, 0.36, 0.48, 0.8, true);
    Rotation reverted = r.revert();
    checkRotation(r.compose(reverted, RotationConvention.VECTOR_OPERATOR), 1, 0, 0, 0);
    checkRotation(reverted.compose(r, RotationConvention.VECTOR_OPERATOR), 1, 0, 0, 0);
    Assert.assertEquals(r.getAngle(), reverted.getAngle(), 1.0e-12);
    Assert.assertEquals(-1,
                        Vector3D.dotProduct(r.getAxis(RotationConvention.VECTOR_OPERATOR),
                                           reverted.getAxis(RotationConvention.VECTOR_OPERATOR)),
                        1.0e-12);
  }

  @Test
  public void testRevertFrameTransform() {
    Rotation r = new Rotation(0.001, 0.36, 0.48, 0.8, true);
    Rotation reverted = r.revert();
    checkRotation(r.compose(reverted, RotationConvention.FRAME_TRANSFORM), 1, 0, 0, 0);
    checkRotation(reverted.compose(r, RotationConvention.FRAME_TRANSFORM), 1, 0, 0, 0);
    Assert.assertEquals(r.getAngle(), reverted.getAngle(), 1.0e-12);
    Assert.assertEquals(-1,
                        Vector3D.dotProduct(r.getAxis(RotationConvention.FRAME_TRANSFORM),
                                           reverted.getAxis(RotationConvention.FRAME_TRANSFORM)),
                        1.0e-12);
  }

[1061, 'src/test/java', 'org.apache.commons.math3.linear', 'RealVectorAbstractTest', 'testWalkInDefaultOrderChangingVisitor2', 1666, 1714]

[1061, 'src/test/java', 'org.apache.commons.math3.linear', 'RealVectorAbstractTest', 'testWalkInOptimizedOrderChangingVisitor2', 1797, 1845]

    @Test
    public void testWalkInDefaultOrderChangingVisitor2() {
        final RealVector v = create(new double[5]);
        final RealVectorChangingVisitor visitor;
        visitor = new RealVectorChangingVisitor() {

            public double visit(int index, double value) {
                return 0.0;
            }

            public void start(int dimension, int start, int end) {
                // Do nothing
            }

            public double end() {
                return 0.0;
            }
        };
        try {
            v.walkInDefaultOrder(visitor, -1, 4);
            Assert.fail();
        } catch (OutOfRangeException e) {
            // Expected behavior
        }
        try {
            v.walkInDefaultOrder(visitor, 5, 4);
            Assert.fail();
        } catch (OutOfRangeException e) {
            // Expected behavior
        }
        try {
            v.walkInDefaultOrder(visitor, 0, -1);
            Assert.fail();
        } catch (OutOfRangeException e) {
            // Expected behavior
        }
        try {
            v.walkInDefaultOrder(visitor, 0, 5);
            Assert.fail();
        } catch (OutOfRangeException e) {
            // Expected behavior
        }
        try {
            v.walkInDefaultOrder(visitor, 4, 0);
            Assert.fail();
        } catch (NumberIsTooSmallException e) {
            // Expected behavior
        }
    }

    @Test
    public void testWalkInOptimizedOrderChangingVisitor2() {
        final RealVector v = create(new double[5]);
        final RealVectorChangingVisitor visitor;
        visitor = new RealVectorChangingVisitor() {

            public double visit(int index, double value) {
                return 0.0;
            }

            public void start(int dimension, int start, int end) {
                // Do nothing
            }

            public double end() {
                return 0.0;
            }
        };
        try {
            v.walkInOptimizedOrder(visitor, -1, 4);
            Assert.fail();
        } catch (OutOfRangeException e) {
            // Expected behavior
        }
        try {
            v.walkInOptimizedOrder(visitor, 5, 4);
            Assert.fail();
        } catch (OutOfRangeException e) {
            // Expected behavior
        }
        try {
            v.walkInOptimizedOrder(visitor, 0, -1);
            Assert.fail();
        } catch (OutOfRangeException e) {
            // Expected behavior
        }
        try {
            v.walkInOptimizedOrder(visitor, 0, 5);
            Assert.fail();
        } catch (OutOfRangeException e) {
            // Expected behavior
        }
        try {
            v.walkInOptimizedOrder(visitor, 4, 0);
            Assert.fail();
        } catch (NumberIsTooSmallException e) {
            // Expected behavior
        }
    }

[1066, 'src/test/java', 'org.apache.commons.math3.ode.nonstiff', 'DormandPrince853IntegratorTest', 'testTooLargeFirstStep', 190, 216]

[1066, 'src/test/java', 'org.apache.commons.math3.ode.nonstiff', 'GraggBulirschStoerIntegratorTest', 'testTooLargeFirstStep', 273, 299]

  @Test
  public void testTooLargeFirstStep()
      throws DimensionMismatchException, NumberIsTooSmallException,
             MaxCountExceededException, NoBracketingException {

      AdaptiveStepsizeIntegrator integ =
              new DormandPrince853Integrator(0, Double.POSITIVE_INFINITY, Double.NaN, Double.NaN);
      final double start = 0.0;
      final double end   = 0.001;
      FirstOrderDifferentialEquations equations = new FirstOrderDifferentialEquations() {

          public int getDimension() {
              return 1;
          }

          public void computeDerivatives(double t, double[] y, double[] yDot) {
              Assert.assertTrue(t >= FastMath.nextAfter(start, Double.NEGATIVE_INFINITY));
              Assert.assertTrue(t <= FastMath.nextAfter(end,   Double.POSITIVE_INFINITY));
              yDot[0] = -100.0 * y[0];
          }

      };

      integ.setStepSizeControl(0, 1.0, 1.0e-6, 1.0e-8);
      integ.integrate(equations, start, new double[] { 1.0 }, end, new double[1]);

  }

  @Test
  public void testTooLargeFirstStep()
      throws DimensionMismatchException, NumberIsTooSmallException,
             MaxCountExceededException, NoBracketingException {

      AdaptiveStepsizeIntegrator integ =
              new GraggBulirschStoerIntegrator(0, Double.POSITIVE_INFINITY, Double.NaN, Double.NaN);
      final double start = 0.0;
      final double end   = 0.001;
      FirstOrderDifferentialEquations equations = new FirstOrderDifferentialEquations() {

          public int getDimension() {
              return 1;
          }

          public void computeDerivatives(double t, double[] y, double[] yDot) {
              Assert.assertTrue(t >= FastMath.nextAfter(start, Double.NEGATIVE_INFINITY));
              Assert.assertTrue(t <= FastMath.nextAfter(end,   Double.POSITIVE_INFINITY));
              yDot[0] = -100.0 * y[0];
          }

      };

      integ.setStepSizeControl(0, 1.0, 1.0e-6, 1.0e-8);
      integ.integrate(equations, start, new double[] { 1.0 }, end, new double[1]);

  }

[1101, 'src/test/java', 'org.apache.commons.math3.linear', 'BlockFieldMatrixTest', 'testDimensions', 160, 170]

[1101, 'src/test/java', 'org.apache.commons.math3.linear', 'FieldMatrixImplTest', 'testDimensions', 102, 112]

    @Test
    public void testDimensions() {
        BlockFieldMatrix<Fraction> m = new BlockFieldMatrix<Fraction>(testData);
        BlockFieldMatrix<Fraction> m2 = new BlockFieldMatrix<Fraction>(testData2);
        Assert.assertEquals("testData row dimension",3,m.getRowDimension());
        Assert.assertEquals("testData column dimension",3,m.getColumnDimension());
        Assert.assertTrue("testData is square",m.isSquare());
        Assert.assertEquals("testData2 row dimension",m2.getRowDimension(),2);
        Assert.assertEquals("testData2 column dimension",m2.getColumnDimension(),3);
        Assert.assertTrue("testData2 is not square",!m2.isSquare());
    }

    @Test
    public void testDimensions() {
        Array2DRowFieldMatrix<Fraction> m = new Array2DRowFieldMatrix<Fraction>(testData);
        Array2DRowFieldMatrix<Fraction> m2 = new Array2DRowFieldMatrix<Fraction>(testData2);
        Assert.assertEquals("testData row dimension",3,m.getRowDimension());
        Assert.assertEquals("testData column dimension",3,m.getColumnDimension());
        Assert.assertTrue("testData is square",m.isSquare());
        Assert.assertEquals("testData2 row dimension",m2.getRowDimension(),2);
        Assert.assertEquals("testData2 column dimension",m2.getColumnDimension(),3);
        Assert.assertTrue("testData2 is not square",!m2.isSquare());
    }

[1102, 'src/test/java', 'org.apache.commons.math3.linear', 'BlockFieldMatrixTest', 'testAdd', 185, 196]

[1102, 'src/test/java', 'org.apache.commons.math3.linear', 'FieldMatrixImplTest', 'testAdd', 126, 137]

    @Test
    public void testAdd() {
        BlockFieldMatrix<Fraction> m = new BlockFieldMatrix<Fraction>(testData);
        BlockFieldMatrix<Fraction> mInv = new BlockFieldMatrix<Fraction>(testDataInv);
        FieldMatrix<Fraction> mPlusMInv = m.add(mInv);
        Fraction[][] sumEntries = mPlusMInv.getData();
        for (int row = 0; row < m.getRowDimension(); row++) {
            for (int col = 0; col < m.getColumnDimension(); col++) {
                Assert.assertEquals(testDataPlusInv[row][col],sumEntries[row][col]);
            }
        }
    }

    @Test
    public void testAdd() {
        Array2DRowFieldMatrix<Fraction> m = new Array2DRowFieldMatrix<Fraction>(testData);
        Array2DRowFieldMatrix<Fraction> mInv = new Array2DRowFieldMatrix<Fraction>(testDataInv);
        FieldMatrix<Fraction> mPlusMInv = m.add(mInv);
        Fraction[][] sumEntries = mPlusMInv.getData();
        for (int row = 0; row < m.getRowDimension(); row++) {
            for (int col = 0; col < m.getColumnDimension(); col++) {
                Assert.assertEquals(testDataPlusInv[row][col],sumEntries[row][col]);
            }
        }
    }

[1124, 'src/test/java', 'org.apache.commons.math3.analysis.integration', 'IterativeLegendreGaussIntegratorTest', 'testIssue464', 124, 158]

[1124, 'src/test/java', 'org.apache.commons.math3.analysis.integration', 'LegendreGaussIntegratorTest', 'testIssue464', 105, 138]

    @Test
    public void testIssue464() {
        final double value = 0.2;
        UnivariateFunction f = new UnivariateFunction() {
            public double value(double x) {
                return (x >= 0 && x <= 5) ? value : 0.0;
            }
        };
        IterativeLegendreGaussIntegrator gauss
            = new IterativeLegendreGaussIntegrator(5, 3, 100);

        // due to the discontinuity, integration implies *many* calls
        double maxX = 0.32462367623786328;
        Assert.assertEquals(maxX * value, gauss.integrate(Integer.MAX_VALUE, f, -10, maxX), 1.0e-7);
        Assert.assertTrue(gauss.getEvaluations() > 37000000);
        Assert.assertTrue(gauss.getIterations() < 30);

        // setting up limits prevents such large number of calls
        try {
            gauss.integrate(1000, f, -10, maxX);
            Assert.fail("expected TooManyEvaluationsException");
        } catch (TooManyEvaluationsException tmee) {
            // expected
            Assert.assertEquals(1000, tmee.getMax());
        }

        // integrating on the two sides should be simpler
        double sum1 = gauss.integrate(1000, f, -10, 0);
        int eval1   = gauss.getEvaluations();
        double sum2 = gauss.integrate(1000, f, 0, maxX);
        int eval2   = gauss.getEvaluations();
        Assert.assertEquals(maxX * value, sum1 + sum2, 1.0e-7);
        Assert.assertTrue(eval1 + eval2 < 200);

    }

    @Test
    public void testIssue464() {
        final double value = 0.2;
        UnivariateFunction f = new UnivariateFunction() {
            public double value(double x) {
                return (x >= 0 && x <= 5) ? value : 0.0;
            }
        };
        LegendreGaussIntegrator gauss = new LegendreGaussIntegrator(5, 3, 100);

        // due to the discontinuity, integration implies *many* calls
        double maxX = 0.32462367623786328;
        Assert.assertEquals(maxX * value, gauss.integrate(Integer.MAX_VALUE, f, -10, maxX), 1.0e-7);
        Assert.assertTrue(gauss.getEvaluations() > 37000000);
        Assert.assertTrue(gauss.getIterations() < 30);

        // setting up limits prevents such large number of calls
        try {
            gauss.integrate(1000, f, -10, maxX);
            Assert.fail("expected TooManyEvaluationsException");
        } catch (TooManyEvaluationsException tmee) {
            // expected
            Assert.assertEquals(1000, tmee.getMax());
        }

        // integrating on the two sides should be simpler
        double sum1 = gauss.integrate(1000, f, -10, 0);
        int eval1   = gauss.getEvaluations();
        double sum2 = gauss.integrate(1000, f, 0, maxX);
        int eval2   = gauss.getEvaluations();
        Assert.assertEquals(maxX * value, sum1 + sum2, 1.0e-7);
        Assert.assertTrue(eval1 + eval2 < 200);

    }

[1147, 'src/test/java', 'org.apache.commons.math3.analysis.differentiation', 'SparseGradientTest', 'testPrimitiveAdd', 137, 142]

[1147, 'src/test/java', 'org.apache.commons.math3.analysis.differentiation', 'SparseGradientTest', 'testPrimitiveMultiply', 169, 174]

    @Test
    public void testPrimitiveAdd() {
        checkF0F1(SparseGradient.createVariable(0, 1.0).add(5), 6.0, 1.0, 0.0, 0.0);
        checkF0F1(SparseGradient.createVariable(1, 2.0).add(5), 7.0, 0.0, 1.0, 0.0);
        checkF0F1(SparseGradient.createVariable(2, 3.0).add(5), 8.0, 0.0, 0.0, 1.0);
    }

    @Test
    public void testPrimitiveMultiply() {
        checkF0F1(SparseGradient.createVariable(0, 1.0).multiply(5),  5.0, 5.0, 0.0, 0.0);
        checkF0F1(SparseGradient.createVariable(1, 2.0).multiply(5), 10.0, 0.0, 5.0, 0.0);
        checkF0F1(SparseGradient.createVariable(2, 3.0).multiply(5), 15.0, 0.0, 0.0, 5.0);
    }

[1148, 'src/test/java', 'org.apache.commons.math3.linear', 'QRDecompositionTest', 'testDimensions', 54, 68]

[1148, 'src/test/java', 'org.apache.commons.math3.linear', 'RRQRDecompositionTest', 'testDimensions', 53, 67]

    @Test
    public void testDimensions() {
        checkDimension(MatrixUtils.createRealMatrix(testData3x3NonSingular));

        checkDimension(MatrixUtils.createRealMatrix(testData4x3));

        checkDimension(MatrixUtils.createRealMatrix(testData3x4));

        Random r = new Random(643895747384642l);
        int    p = (5 * BlockRealMatrix.BLOCK_SIZE) / 4;
        int    q = (7 * BlockRealMatrix.BLOCK_SIZE) / 4;
        checkDimension(createTestMatrix(r, p, q));
        checkDimension(createTestMatrix(r, q, p));

    }

    @Test
    public void testDimensions() {
        checkDimension(MatrixUtils.createRealMatrix(testData3x3NonSingular));

        checkDimension(MatrixUtils.createRealMatrix(testData4x3));

        checkDimension(MatrixUtils.createRealMatrix(testData3x4));

        Random r = new Random(643895747384642l);
        int    p = (5 * BlockRealMatrix.BLOCK_SIZE) / 4;
        int    q = (7 * BlockRealMatrix.BLOCK_SIZE) / 4;
        checkDimension(createTestMatrix(r, p, q));
        checkDimension(createTestMatrix(r, q, p));

    }

[1153, 'src/test/java', 'org.apache.commons.math3.random', 'EmpiricalDistributionTest', 'testLoadNullURL', 224, 227]

[1153, 'src/test/java', 'org.apache.commons.math3.random', 'EmpiricalDistributionTest', 'testLoadNullFile', 229, 232]

    @Test(expected=NullArgumentException.class)
    public void testLoadNullURL() throws Exception {
        new EmpiricalDistribution().load((URL) null);
    }

    @Test(expected=NullArgumentException.class)
    public void testLoadNullFile() throws Exception {
        new EmpiricalDistribution().load((File) null);
    }

[1165, 'src/test/java', 'org.apache.commons.math3.transform', 'RealTransformerAbstractTest', 'testTransformReal', 253, 263]

[1165, 'src/test/java', 'org.apache.commons.math3.transform', 'RealTransformerAbstractTest', 'testTransformFunction', 276, 286]

    @Test
    public void testTransformReal() {
        final TransformType[] type = TransformType.values();
        for (int i = 0; i < getNumberOfValidDataSizes(); i++) {
            final int n = getValidDataSize(i);
            final double tol = getRelativeTolerance(i);
            for (int j = 0; j < type.length; j++) {
                doTestTransformReal(n, tol, type[j]);
            }
        }
    }

    @Test
    public void testTransformFunction() {
        final TransformType[] type = TransformType.values();
        for (int i = 0; i < getNumberOfValidDataSizes(); i++) {
            final int n = getValidDataSize(i);
            final double tol = getRelativeTolerance(i);
            for (int j = 0; j < type.length; j++) {
                doTestTransformFunction(n, tol, type[j]);
            }
        }
    }

[1170, 'src/test/java', 'org.apache.commons.math3.linear', 'BlockFieldMatrixTest', 'testGetSetRowLarge', 1040, 1059]

[1170, 'src/test/java', 'org.apache.commons.math3.linear', 'BlockFieldMatrixTest', 'testGetSetColumnLarge', 1103, 1122]

    @Test
    public void testGetSetRowLarge() {
        int n = 3 * BlockFieldMatrix.BLOCK_SIZE;
        FieldMatrix<Fraction> m = new BlockFieldMatrix<Fraction>(FractionField.getInstance(), n, n);
        Fraction[] sub = new Fraction[n];
        Arrays.fill(sub, new Fraction(1));

        m.setRow(2, sub);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                if (i != 2) {
                    Assert.assertEquals(new Fraction(0), m.getEntry(i, j));
                } else {
                    Assert.assertEquals(new Fraction(1), m.getEntry(i, j));
                }
            }
        }
        checkArrays(sub, m.getRow(2));

    }

    @Test
    public void testGetSetColumnLarge() {
        int n = 3 * BlockFieldMatrix.BLOCK_SIZE;
        FieldMatrix<Fraction> m = new BlockFieldMatrix<Fraction>(FractionField.getInstance(), n, n);
        Fraction[] sub = new Fraction[n];
        Arrays.fill(sub, new Fraction(1));

        m.setColumn(2, sub);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                if (j != 2) {
                    Assert.assertEquals(new Fraction(0), m.getEntry(i, j));
                } else {
                    Assert.assertEquals(new Fraction(1), m.getEntry(i, j));
                }
            }
        }
        checkArrays(sub, m.getColumn(2));

    }

[1172, 'src/test/java', 'org.apache.commons.math3.random', 'Well19937aTest', 'testReferenceCode', 29, 113]

[1172, 'src/test/java', 'org.apache.commons.math3.random', 'Well19937cTest', 'testReferenceCode', 29, 113]

    @Test
    public void testReferenceCode() {
        int[] base = {
              740849862,  1202665156,  -199039369,  -259008301,  -291878969, -1164428990, -1565918811,   491009864,
            -1883086670,  1383450241,  1244617256,   689006653, -1576746370, -1307940314,  1421489086,  1742094000,
             -595495729,  1047766204,  1875773301, -1637793284,  1379017098,   262792705,   191880010,  -251000180,
            -1753047622,  -972355720,    90626881,  1644693418,  1503365577,   439653419,  1806361562,  1268823869
         };
        int[] init = new int[624];
        for (int i = 0; i < init.length; ++i) {
            init[i] = base[i % base.length] + i;
        }
        Well19937a mt = new Well19937a(init);
        int[] refInt = {
            -612874471,   -354976292,  -1838197125,  -1781560577,    278390997,   1214938280,  -1752615390,   -760835246,  -1712883765,   -241205782,   -145390202,    495649160,   -514388259,  -1271015916,  -1640000013,    849273623,
            -549729394,  -1206917255,   -545909692,    811925434,  -1665729633,  -1525292882,   1416246482,   -153220826,   1148868872,   -326143196,   1724979062,   1790931148,  -1648679618,   -439051683,    112482777,  -1484051520,
           -1881272572,  -1270447031,  -1216102401,   1579107248,  -1224395621,   2144411988,   -416216641,  -1529222361,   1628987080,    164445245,   1506928916,    928145916,   1436000427,    862025970,    560077705,  -1887251027,
            -514360858,   1735094506,    475624879,   1755802355,    295448361,   -155399225,      3972415,   1368201076,   -465126094,  -1622687259,   -246099304,   1798631152,  -1937269102,  -1700560396,   -293352622,   -896632303,
           -2088933220,   -194382452,   -480218162,  -1618517785,  -1925031481,   -150217434,   1678937261,   2130832364,   -485546678,  -1499224981,   1430390884,  -1895417302,    210514746,   1781140999,  -1940853105,  -1238099647,
             485922557,   -103223212,    633481679,   -632946979,    695235541,  -1661735272,    277603567,   -958341538,    256982285,   1850270018,   -327388076,   -219053874,   1380560653,  -1221689980,   1335863752,   -545032383,
            -575291735,  -1295623907,   -140058298,   1063302709,  -1290702617,   -790401546,   -170630961,  -1203114473,   1458063108,  -1212753301,   1546428514,   2112636413,  -1463028928,  -1812598032,   -883529486,   1131084094,
              62042165,   2135819802,  -1192342739,     98361522,  -1341042205,   -475283063,  -1632033747,   1745196892,    168608689,   -914987039,    274428907,   -881357258,    167940012,  -1975737532,   -903960486,  -1370984244,
            -589352935,   1783633514,   -570111010,     71495377,    194463285,  -1243905021,  -1398490898,    221691209,    -55728834,   -638916786,   -770622372,  -1911651810,   -295233027,    301467998,   2058638784,    681490183,
           -1547865078,  -1668135684,   1299261826,   1649468635,    287995017,  -2076844852,   1193468826,   -853948258,    120082777,   1051829542,  -1288514343,   -159456430,    275748820,   -480127107,   -604943233,  -2138088332,
            1202614819,   1427201263,  -1906344469,  -1230779533,   1690367192,    733159097,    794410312,  -1114452505,  -1601554413,    976747949,   1517787154,   2091780205,   1052078906,   1919282771,   -191013374,   1805397142,
             736939268,  -1056272823,   -727464316,   -659459005,    797803875,  -1104633884,   1042342081,    -24514837,   1919469940,   1903722546,   -814157872,   1605407665,   -262351256,   -288949635,    729204844,  -1132605534,
             745453338,    387915035,   1094173337,   2100279147,    156863702,   -257377544,   -719587984,  -1496015613,   1908993744,   2016957554,    918749666,   -135963651,  -1356808639,  -1711185741,   1472589240,   -398100149,
             628791415,  -1381837652,  -1820702771,   -593586943,  -1456631279,  -1837975351,  -1394249972,   -556916726,    833231177,     43449750,   1029237092,  -2086437337,   -459463076,   -533031784,  -1739648287,  -1374722961,
            2024908394,   1389678488,      2018558,  -1391707864,   -795935743,    904816957,    836583280,   1766194531,  -1374431014,   -904437876,   2030248636,   -265724199,   2056758426,   -810499837,    887193593,    -77811488,
            1496312336,  -1874348275,   -456193866,  -2137130942,    868120387,     29025455,  -1999867716,   2001322335,   -579152815,   -390892056,   1592011837,   -306394879,     93636886,   -190879994,   1923358153,    269052141,
            -396050253,   -987531729,    480350991,   1276744541,  -1445571957,   -957571005,  -2046270221,  -1715395752,   1113585628,  -1782113514,   -697560146,    835320000,   1014320959,  -2119834109,    460056841,  -1464772991,
           -1282790418,  -2120806165,     86176097,   -731086307,    832497517,  -1876684928,    541008240,    551124479,   -450919132,    647860281,  -2115397586,    979247589,   1095559204,   1927958688,    169497703,   1999579054,
            2019745038,   1656022059,  -1109662138,    375237154,   1450814436,    919988416,    849761266,   1457057327,   1771166577,  -1639880487,   -852488298,   1767063646,    657295386,   -585561879,    740792583,   1664558308,
            -654749506,   1109275990,    182597559,   1106789745,  -1806628480,     25948116,   1748374299,    196057325,   -164213209,   1687024594,    782029276,   1879737947,  -1528219611,    412585737,   1190325629,   1985821911,
           -1272945202,  -1238637137,    465818730,  -1537670961,   1131953615,    905623579,    609183424,   1138422991,   1522974699,    589719061,  -1310894604,    890952933,   -885204790,   -393535694,   1238408670,   1780660354,
             677803525,  -1121509064,   1553148616,   1109165936,  -1450120385,   1525252521,  -1354897489,   -595402189,  -1274551767,   -869281409,   1788815975,   2020262116,   1124100185,   -400839020,    310574108,   1354413045,
           -1310514485,   1895732085,    626639054,   1667355357,   2065637178,  -1889009143,   -440157749,   1762849463,  -1693853642,    -56602956,   -930874188,   -398470740,    778356402,  -2113156881,     42854964,   1844399604,
           -2098310302,  -1812029757,   1441188713,    899579267,   1266994172,   1841370863,   -660740252,    -43254718,   1124500192,   1884907320,    879997211,   1775139845,  -1360112721,   1630490057,    362567879,   1113475029,
             290319279,  -1209506867,    398146039,   -957742350,   1185761854,   1519676447,   -912689915,  -1117128973,   -305563462,  -1928033363,  -1766324543,   1702753492,   1696951912,  -1895072395,    932663591,   -566548128,
             991675996,     56529814,    980735023,    718166662,   -650028466,   -886842051,   1857048587,   -569023569,  -1820572202,   -851452711,   -958700452,   -621825633,    -65649888,   -510143183,    761267599,  -1692108035,
            1729071710,   1623630864,    -53498654,    267235687,    659201413,   1152882627,   -824194574,    356636960,   -502391121,   -538453360,     66115376,  -1633290370,  -1522088932,    268949070,    684499443,   -859474501,
            1586764345,  -1515639709,    319695602,   -307025150,     69076508,   1050726785,  -1340930110,    552191600,   -207852941,   -273572993,   -539580440,    710343120,   1957076127,  -1107172811,   -561518280,  -1775022699,
            1978792904,   1935531451,  -2084046304,   -419742902,   -737652926,    614022023,   1676952428,    769892939,  -1092786807,  -1117113223,   -266029995,   -350150999,    207738542,   1964896575,     48805284,   1736500159,
             551289617,  -1847923501,   1856609505,   2007480480,   -681860219,  -1198106493,   1483591043,   -523895316,  -1814473078,  -1521087404,  -1348859926,   1298056897,   1813789478,    946683654,     79196400,   1766250931,
             472737685,   1764634332,  -1844726079,   -130619045,   -508713868,  -1762537125,   1010108863,    170107098,   1705386380,  -1139681802,    183739097,   1662699401,   1842694501,   1714633805,     46208876,    616720693,
            -252553427,   1986302230,   -103130254,   1943225981,    110746655,    553260552,   1588938073,  -1934623163,  -2144781332,  -2086217416,   1941265852,   -781953226,   1216234254,    605543697,   -710872598,   2048636577,
           -1986927728,  -1007017623,   1243051501,   -614249563,  -2128221291,    581579813,   1173464240,  -1906830937,    261329601,  -1805974103,    769823490,   1858731164,   -561762071,    516417430,  -1221329437,   -825500715,
            1091364656,   -993658663,  -1475434188,  -1070804384,  -1876492082,    899494424,    683486936,    878807455,     56642807,  -1268202879,   1379172046,  -1386869373,  -1158233876,   1759190552,   1597629789,   1411151497,
           -1254268471,   1075936979,   -918778269,  -2132675184,    953140888,   1906982077,   1154200766,   -365384600,  -1142488826,    708535121,  -2134869964,  -1531201665,  -2100526761,   1268256467,   2071480803,    193135243,
            1374158182,    989505347,   -933612202,  -2134839213,  -1302795271,  -2092029041,   1812014826,   2090855917,   2005348528,    606434393,    -60141386,     11156360,    539516285,   -122485034,   -893237911,   -978127424,
            1509901816,   -451029719,    428544700,  -1622965963,  -1993611605,  -1989324583,   1104111587,   -795138585,   -899552401,  -2110167769,   -234502445,   1586963605,   -503778455,    529261062,    325327284,   -106186403,
              65369563,  -1475700698,   -228624261,    715975009,   1099352363,  -1796883396,   1376542700,   -308942420,   -344940451,   -395389249,  -1562737166,   1869802677,   1273494710,   2075587668,   -789570273,   1563347596,
            1142901755,   1676422422,  -1729157809,  -1399423717,  -1814262429,  -1809707284,   1393992342,   -570246212,   1065528749,   -781643849,   1218667301,  -1097949471,   1305629790,    901301039,   -704762030,    360582612,
            1411910672,   1848068741,   -614500891,   -146889637,   -913903597,    723527277,   -147033328,   -199273155,    734997691,  -2072735286,   2129258691,  -1385074104,    931616624,   1065477319,  -1543474555,   -531410292,
           -2123119121,  -1538464113,  -1153585193,   1559931968,   -654877011,    879865200,   1489681397,   1998864644,  -1964160144,    163671782,   -858364148,   -323324233,    801208648,    705864113,    436184243,    643773864,
            2087594507,    134637265,   -749956494,  -1657343972,  -1828172168,    -27357303,  -1145161336,  -1192513644,    216148260,    611393153,    -13752671,   -358631090,  -1211920749,    593572064,    657629904,  -1445961088,
            -250704995,   1797542707,  -2122311891,   -316774825,   -296303057,   -868002056,    -86697533,   2020588145,   1203427903,  -1371839056,    669531557,  -2031033836,   1323994690,     13703036,    785437772,  -1465821554,
            -694756014,  -2131068154,  -1745448876,  -1095891733,    936594025,  -1119068454,    855423970,   1705079340,   -905640608,    162297141,   1336619311,   -344353769,    -92608588,  -1080573824,   2002293105,  -2088030765,
           -1684198727,   -129054718,   -949437132,   -127983221,   -216664110,   1700146143,   -711174649,   1500113839,   1212236226,  -2017364219,  -1263597675,    511929344,   -323998524,  -2021313185,   1803000924,    927670608,
             336267187,   1244256964,  -1665390108,    991395134,   -251232188,   1267445783,   1547951569,    740269916,   1776431169,   1687220659,    228229817,    271386089,   -682906779,   -438090344,   1190436796,   -744272540,
            1879221151,   1145200306,  -1730983338,  -1119209309,     90826726,   1567861540,   1638852830,  -1645384932,   1566909531,   1088584561,   1030555565,  -1872052014,    720320695,   -885053674,   -321216789,    739907579,
             368580703,   -443635520,   1232705619,  -1355949988,  -1047211249,  -1571429448,    599299852,   1036970439,   1513838571,    -51797291,    -26647565,  -1262878942,   -916262582,   1579082269,   -292007383,   1289013866,
           -1612184284,   1451738668,    448608569,    476432992,  -1229609565,    786372409,    929928149,   -150100614,    448155944,  -1322320576,   -856549627,   1057443268,  -1536809554,    577508258,    584906122,    275295163,
            -604262071,   -236043234,  -1866434954,  -2072447013,    646132876,    847562546,   -310005953,  -1104162658,    393261203,   -730102354,    440824482,   1654535035,  -1296359745,   1487359328,   -977776604,   -775827779,
           -1298695106,    519080622,   1697162240,    227873031,   -371123123,   1273702312,  -1710063656,  -2138342344,   1139555478,   1531578907,  -1498880699,   1183507362,   1875307493,  -1649740413,   2135386504,   -962458407,
             424161768,    504272962,    202204247,   1783466420,   2015579232,   -676642965,   2067456450,    914480415,   -620398841,   1880405399,   1406637142,   1951104977,    633496157,    224861869,    -58659291,    994942775,
            -479000645,   1421449115,    100168104,    249754169,  -1219011494,   1736303638,    364013694,  -1750035055,   -479217141,   1652913106,  -2109452331,   1633842910,  -1547663337,    936627493,  -1152799743,    896955899,
           -1407742850,   -523769014,    357161414,    872293304,    744895980,    720829676,   -240843156,   -111779524,   1292836315,  -1792141538,   1946959925,   1181751089,  -1120674052,   1185192575,  -1387002557,   1973209255,
            -120887476,   -766577735,   -443913073,    786620227,    428564781,   -101232106,   -425959852,    198082021,   1173272226,  -1744840378,  -1621135606,  -1539498583,  -1101274572,     43399711,  -1256764602,   1201920787,
            2049426139,    846545551,  -2121520873,  -1202939675,   -470425740,    321987390,   1862019060,   -951540342,   -894238318,   -430407175,  -1662746491,    656574776,   1580373777,   -431290218,   1645824323,  -1953526979,
            -374682356,    474291752,   1071558425,    511038981,   -760598678,   -567797285,  -1176476266,   -268409005,  -2130644484,    -67970563,   1756046948,   1429860462,  -1130984739,   -124916495,  -1544436836,  -1863524031,
            1024916487,  -1388636482,  -1573205065,    892628956,   1831270021,   1176430590,   1158914682,  -2006787098,  -1228130033,   1516111488,  -1499151347,    470546266,   1642603981,   1425140838,  -1823071475,  -1775267236,
           -1009380612,    164746986,   1129677098,   1842642579,   -482342932,   -507480364,   1656012309,   1981601761,   -881042120,   -511987083,    342447017,    381192578,    983008095,    741012865,  -1877136350,   -199211983,
            -452784912,   1929572576,  -1678291139,   -864375281,  -1610561247,  -1936356726,   -749553767,   -865893512,   -567081879,  -1303973729,   -939636958,   -622974563,    428284937,   1049237414,    852280765,     86648946,
           -1353851401,  -1045422335,    898035731,  -1636093996,  -1083174191,    245046915,   -359768226,  -1028491655,   1051575118,   1774289451,   1839389415,  -1594053468,    736707953,   1873556950,    401186168,   -583669552,
             -88375334,   2002752071,    264506453,  -1304812107,   -759203942,   -114958524,  -1878903503,    841613720,   1910863820,  -1738114003,    701455920,   1791058048,  -1850960547,   1672292671,   1172188809,    604848896,
           -1607489375,    305370478,   -948153885,  -1971080100,  -1848966954,   -584538365,     39416319,  -1689119162,    944942598,   1777111075,   1534005553,   2022718432,    -25820385,      3077695,   -315950520,   1859184648,
           -1397829266,  -1666371809,    858913807,   -610818620,   1554973298,    580023809,  -1662988256,   -408630026,   1316681876,    738204271,    942829881,   -758486983,    780345857,    667165037,  -2086803585,    789741324
        };

        for (int i = 0; i < refInt.length; ++i) {
            Assert.assertEquals(refInt[i], mt.nextInt());
        }

    }

    @Test
    public void testReferenceCode() {
        int[] base = {
              740849862,  1202665156,  -199039369,  -259008301,  -291878969, -1164428990, -1565918811,   491009864,
            -1883086670,  1383450241,  1244617256,   689006653, -1576746370, -1307940314,  1421489086,  1742094000,
             -595495729,  1047766204,  1875773301, -1637793284,  1379017098,   262792705,   191880010,  -251000180,
            -1753047622,  -972355720,    90626881,  1644693418,  1503365577,   439653419,  1806361562,  1268823869
         };
        int[] init = new int[624];
        for (int i = 0; i < init.length; ++i) {
            init[i] = base[i % base.length] + i;
        }
        Well19937c mt = new Well19937c(init);
        int[] refInt = {
            2128528153,    327121884,    935445371,    -83026433,  -1041143083,   2084595880,  -1073535198,  -1678863790,   -523636021,  -1514837782,   -736786810,   1527711112,  -1051227939,    978703380,    410322163,   1727815703,
            -648426354,    636056441,   1954420292,     17754810,   -958628705,  -1091307602,   1793078738,  -1680336346,   1792171272,    941973796,  -2066152330,  -1248758068,  -1061211586,    262020189,   1276960217,   -233886784,
            1767509252,  -1811939255,   -406116097,   -742435920,  -1349799525,    240329556,   -332161345,   1488943143,   -332244280,   2093328957,    674753300,  -1930135556,    257111467,     63793650,  -1964335223,   1315849133,
            -797349146,   1372022250,  -1451892049,  -1325138957,   -870401239,  -1294317369,     91490879,    386205044,   -704074702,  -1230679067,   1513674392,   -262996240,   1196007314,   1398903796,    803719762,  -1750926831,
           -1268814180,   1233515404,   1498313934,   -970591257,    611113671,   -261632474,   1834097325,   1709440492,   -150396854,   2120561003,    -62645660,    479080234,   1535125050,   1823378695,  -1129289329,  -1095198399,
            2092564733,     78836308,   -692015409,   1647147229,  -1847922219,   1838279320,   -848333841,  -1375151778,    920238861,   1512628290,   -749439404,    288851918,   -427218675,    679640964,    425700808,  -2077624511,
           -1929434455,   -647176419,    650437190,  -1926749131,  -1564744729,    734494454,    108193743,    246246679,    810042628,   1952337771,   1089253730,  -1874275331,   1428419392,   -492969232,   1945270770,   -201265602,
            -755490251,   -624426214,   -699605715,   -113446478,    809091299,  -1521531511,   1136505389,   -523660964,    132928433,   1926559713,  -1485314325,   -508322506,     46307756,  -1627479740,   -589386406,  -1855555892,
             584299545,   1272841066,   -597242658,    925134545,   1102566453,   -753335037,     -9523218,  -1778632375,    568963646,    764338254,   1259944540,  -2000124642,   1307414525,   -151384482,    807294400,   1993749511,
             -15503094,   -709471492,   2104830082,   1387684315,  -1929056119,    224254668,   -733550950,   -889466978,  -1987783335,   -437144026,    995905753,  -1021386158,  -1096313388,  -1014152835,  -1303258241,   1201884788,
           -1845042397,   1421462511,    980805867,   2143771251,    481226968,   1790544569,    328448328,   1995857639,    -66668269,  -1411421267,   -222586606,    866950765,   -308713926,  -1048350893,    993222402,  -1139265642,
            -871837948,   1145571913,    381928580,     35386691,   1640961123,  -1192981020,    775971009,    594246635,   1603197812,   -575106766,   2023682000,  -1636301903,   -718093720,  -1666421635,  -2146115988,    320593570,
             287355418,    454400027,   1112753817,   1751196267,    782077910,  -1478447368,  -1007557264,   -862315517,  -2035355952,   2123515250,   -557641502,  -1789932035,    879640129,     44167603,    791148984,   1382939723,
           -2135684233,   1825489580,    937345485,  -1839983359,  -1536880111,  -1472578359,   1548052748,  -1471535862,    -14508727,   1509621398,  -2134967452,   -787485401,    815341660,   -327905128,   1028096737,    866906991,
           -1585990806,    859229080,    234806270,    998518056,  -1897890815,   -900923587,   1179856752,   1529572451,    620486106,   1119836556,   1661285564,   2097404633,  -1437490790,    265306115,   -984880135,   1326751968,
            1280043536,    680210701,    155786166,   1550973250,   -325781949,   -597789777,     -1939780,   1345275487,   1930450001,    941449704,    669301309,    693651713,   -990721514,    582968326,    976132553,  -1892942099,
           -1065070157,   -711990993,   -688974833,  -1026091683,   1115346827,  -1305730749,  -1733626381,   -364566696,    -21761572,    -37152746,   -262011730,   1302722752,  -1806313409,   -767072509,    764112137,   1671157377,
            1837645038,  -1021606421,  -1781898911,   -232127459,   -310742675,  -1818095744,  -1128320656,   -705565953,   -354445532,   -523172807,   -433877202,    131904485,    -64292316,    381829280,    229820263,   1797992622,
            1359665678,    978481451,   -885267130,  -1415988446,   -356533788,   -961419072,   1938703090,    708344111,    679299953,    744615129,   1328811158,   1257588574,    569216282,   -753296151,  -1519838713,   2016884452,
            1062684606,   1561736790,   2028643511,  -1353001615,    886376832,   1466953172,   1664783899,   1290079981,    -57483993,  -1176112430,   1634916316,   1976304475,   1374136869,   -648738039,   1058175869,   -909000745,
           -1526439218,    726626991,   2066596202,     64980943,    -26166577,   -885900005,  -1821546816,  -1103727665,    730606315,  -1324948459,   -696956940,  -1300869403,   1171578314,    797249074,  -1600611618,   1928247682,
             307164165,  -1482476232,  -1886179640,   1306433392,   1945271359,  -1272113751,  -1285984081,  -2057145549,    795047465,   1262569087,  -1239828121,   1426641636,   -786371495,   2120199316,   1273690652,     74457589,
           -1033394229,    338952565,     46122958,   1225741533,   2115308090,    678200841,  -1618264885,   -101162569,  -1628976330,  -1232839500,    468709044,   1876019116,     92723122,    233398255,   -554960844,     38494196,
            -406437278,   2083528643,  -1106878615,   -340722557,  -2123964932,    223183343,    108918116,  -1014629054,   -901344544,   -838896840,  -1908460517,  -1763508731,   -926890833,   1703791049,   -667755577,   1694418389,
             791641263,   1095689677,   1119202039,  -1419111438,  -2012259010,    188017439,  -1775110395,  -1971099661,  -1688113734,    131472813,   -776304959,   1511388884,   2080864872,  -1733824651,   1992147495,   1119828320,
            1065336924,  -1357606762,    462963503,   1180719494,   -202678962,   -892646595,    605869323,   1144255663,    878462678,  -1051371303,    872374876,    631322271,   -172600544,  -1552071375,  -1939570033,    151973117,
            1640861022,    310682640,     34192866,   2057773671,  -2004476027,  -1879238973,    582736114,    900581664,   -427390545,  -1232348528,   -535115984,   1321853054,     69386780,  -1729375922,   1418473715,   1022091451,
             496799289,    -80757405,  -1903543310,  -1128846846,      1703964,   1984450945,    856753858,   -812919184,    775486323,  -1376056193,    638628840,    314243536,   1030626207,    644050997,     73923896,    362270613,
             236584904,   1463240891,   -223614432,    435371594,   -751940030,   -124274553,  -1991092884,   1579624267,   1249632649,    157589625,   -345229739,   -366245207,  -1399995986,   1651729983,   1965074340,  -1108970305,
            1163690769,   1732013523,  -1461252895,    669755552,   -476503675,   -264578685,    -32813949,    288222188,    -25734262,    106040916,   1654395626,   -365148479,   2014455846,  -2040447994,   1351639280,   -919975757,
           -1970412139,    -47306532,    222377665,   -363434917,  -1091717516,   2090685531,  -1221091649,  -1729649190,  -1239406708,   1064945398,   -105437479,   -419675255,     74701669,    -12862899,   -498269844,   1566898997,
           -1872838355,   1596887574,    485902962,    469225597,   -881763553,   1307841032,  -1642872487,   1388543045,    379792876,   1095683384,    840780732,   1934378038,   1851278350,  -1359389423,    130868458,   -313448799,
            -663624816,   1031714153,   -608443411,   -205137499,  -1849464427,   1973593637,   1068741808,  -1420655961,   1188762305,    954044841,   -995454462,  -1818101092,  -1937201943,   -324541290,  -1520603933,    572873173,
            -554764496,   1051557081,  -1245136076,   -985349536,    329320398,   1787901464,    -37803304,  -1759310177,  -1463492617,  -1861729663,   1251768782,    256937091,   -779036948,  -2049893864,   1256022877,   1228075657,
           -1550195255,   -611319853,   1190797155,   2047604112,   -576077160,  -1532843331,  -1324899394,   -159729560,   -622525946,  -1080302767,   -236033484,   1895243903,   -410123689,  -1944154157,   -681781021,   1208453003,
             579595878,   1303914051,   -145607082,   -131567277,  -1917288455,    894217359,   -175688726,  -1585480723,    663691440,  -1140068263,   -641711178,   1596080008,    629197693,    976422358,  -1570451095,    525923776,
             895046136,   -504151767,   1602553020,  -1233054923,  -1798474837,  -1488857895,   1055782627,    261863143,   1879276655,    488240679,   1910982611,  -1919441259,    370435945,   1265230086,  -1293284428,  -1503576227,
            2076963035,  -1379628250,   1157098875,   1984461153,  -1947837397,   1705880124,   1453607404,  -1233649748,   1479943773,   -863878721,   -862415630,   -736723275,    940306358,  -1596000684,  -1174889953,   -615723892,
            -885006597,  -1796723178,   1844159055,   -188942309,   2107251811,  -1675486996,  -1009475178,   -859263556,   -431866963,     -9593673,  -1878920923,   -104853791,  -1535224994,    -69315537,    586690130,  -1292234796,
            1378749456,   -301873019,   -319297563,   1677205851,    292450579,  -1289441171,   1788113680,   1907606333,   1464711611,  -1372023606,  -1978832445,  -1772259768,   1949124464,   1818322887,  -1138036603,   1249727628,
           -1474866449,  -1868013169,  -1384567593,    717007936,    954189997,  -1900561040,    738470389,   -158973180,   1732860784,   1936031206,  -1133354740,  -1173166665,   1432976712,    852636081,   1732064691,  -1831788120,
            1273933579,    455403217,   1988395890,    106493468,    506092152,   -610530423,   1698053512,   1311747476,   1969503012,  -1887461759,   1613543073,    903200334,   -737865837,    325656800,  -1234001200,   1492148864,
            2009861533,   -368262605,   1091338541,   2076108119,   -961392337,   1835877112,    316250307,   -853333391,  -2125443777,    815363504,   -798707803,   -158146540,    690786114,   -530775684,   1203556940,   1611485582,
           -1661412270,    -53184506,   2126287444,   -232222229,   1559486057,    283532250,   1202760418,    932144172,   1082594656,   -570104011,    413509167,   -995027177,   -996477516,      -540544,   -745537167,   -712135469,
            -996294983,   -592787198,   1889840948,   1314628747,   -394266926,   -682316577,    456447239,   1728806063,   -396279614,    -43387643,   1915717013,   -861574144,  -1078710588,   -561401249,   1111464540,     63643984,
           -1693870413,   -968369980,  -1053148188,    708799038,   1883537988,    373371671,   -156410415,  -1596483236,  -1846890431,    888692915,  -1025632583,  -1666477591,   -343066267,  -2059058792,    641501628,  -1744347292,
            1648632991,   1743540146,   2020952406,    164014499,    990508262,   1706408228,  -1236471842,   -347116260,   1843634523,    827255665,    300519853,  -1265974830,   -547247177,   -583064554,  -1995437077,    689210107,
             -93151393,    835365056,   1706367315,  -1605902756,    200954895,    431093688,   -277573364,   -928486713,   -552221973,    145432789,   1128919795,   1675095586,   1930359882,   1215849501,  -1447770583,    657776490,
            1885869860,  -1629237204,   -868897479,  -1258169760,   1828140195,   -883850439,    463933909,   -347361158,   1478116648,    801176896,  -1501915899,   1017335748,  -1654508882,    123994786,   1588785290,    791166651,
           -1523108535,    340411166,   -496474762,  -1189711141,     -7392628,   2045171250,  -1245366209,    834787230,  -1346883181,   2034209454,    737043362,    898803323,   1983089087,  -1845404320,      9585188,  -1180608323,
            1665100606,   1949474222,   -211115008,   1151308295,  -2132174259,    913126312,  -2085061672,   1419864120,  -1134542954,    -53833957,   -246913211,    468382370,  -1759479323,   1136686211,   1307012488,  -2036299559,
           -1346099736,    314743106,  -1683101865,   -947151948,   -234529696,  -2103334293,   -279256894,     -1484257,  -1053953017,   1801205399,    941594454,   -874119215,   -672865187,    762284205,  -1494975451,    486607927,
            -898264389,  -1711861093,   -212572760,   2106484281,  -1610786470,   1352525590,   -837779586,   1568282001,   -593019125,  -1146260782,  -1595879979,   -640781858,   1107692311,   1547132709,  -1928385535,  -2057772805,
             634887038,    329772618,    942136006,   -864405576,    501883884,   1537141484,  -1180626836,   1123055420,   1090885851,    421662750,   2033111605,   1710917425,  -1118058244,     74321279,    257328195,  -1199940697,
             208625996,   -442341447,    808119183,   1166827075,   1177417517,  -1856155370,  -1464837036,    -60624923,  -1306220638,    -91104698,  -1434621430,    548899241,     37351476,   1478278431,  -1255061434,    248470035,
            -104642597,  -1865169521,   1418373655,  -1660810523,  -2129015436,    154612798,    276575732,   1930338442,    179503250,   -929294855,    -39452027,  -1377657544,   1442322193,   1137511318,   -432158653,   -984801987,
             743099148,  -1118893528,   -904123623,  -1273146363,  -1884800406,   -803169061,   1254123158,   -484252077,    317646844,    404246525,  -1230293916,   1121445742,    -19657507,    652967153,  -1055406692,   -468950719,
           -1493532921,  -1447624258,  -1369679689,  -1517000228,   -145853307,   1518006526,   1591195514,  -1475557146,   -909722097,   2103182976,   -406830579,  -2124025254,  -1804819507,  -1357512858,    567321869,    409048156,
             567805180,   1749009386,   1762759722,  -1770324077,   1271140844,    468219092,    955792405,   1911965665,   1876314424,   -718200715,  -1278883927,   1392281730,   -120519585,    851473793,    245054754,    -33369039,
            -284877584,   -479534880,   -212346563,   -122017521,  -1461429983,   1331007370,   1788621721,   1739036536,   1446350953,  -1985448033,    685528610,  -1386434659,   1368233993,   2021786790,   1596478397,  -1716635278,
           -2011083017,    171876097,   -311529197,    687812052,    377000657,  -1055547517,  -1499047842,  -1818434951,   -120863666,     33888043,  -1387509273,   -541540700,   1162597745,  -1331415338,   1931708792,   -850270000,
             663845594,   1536495943,   -322924971,  -1380272203,    261190298,   -204874428,  -2104974031,    883819928,    155808204,  -1454446035,   1323388464,  -1696505728,   1549800285,   1018150463,  -1327715703,  -1582480640,
            1013659809,  -1820360082,   1666498787,   1406120540,   -196541482,   1248470531,  -1250433281,    836375878,    177646854,  -1927020253,   2145878321,    689712096,   -596605921,    348283199,   1916993096,    481356808,
            -339687826,   1219340319,    718895887,  -2007521340,  -1859185806,   2042164737,    -58146784,    742449142,   1930754708,    780832111,    715056441,  -1393886151,     -8150527,   -599607443,   -537300865,  -1212516084
        };

        for (int i = 0; i < refInt.length; ++i) {
            Assert.assertEquals(refInt[i], mt.nextInt());
        }

    }

[1184, 'src/test/java', 'org.apache.commons.math3.ode.nonstiff', 'AdamsBashforthIntegratorTest', 'testMinStep', 56, 74]

[1184, 'src/test/java', 'org.apache.commons.math3.ode.nonstiff', 'AdamsMoultonIntegratorTest', 'testMinStep', 56, 76]

    @Test(expected=NumberIsTooSmallException.class)
    public void testMinStep() throws DimensionMismatchException, NumberIsTooSmallException, MaxCountExceededException, NoBracketingException {

          TestProblem1 pb = new TestProblem1();
          double minStep = 0.1 * (pb.getFinalTime() - pb.getInitialTime());
          double maxStep = pb.getFinalTime() - pb.getInitialTime();
          double[] vecAbsoluteTolerance = { 1.0e-15, 1.0e-16 };
          double[] vecRelativeTolerance = { 1.0e-15, 1.0e-16 };

          FirstOrderIntegrator integ = new AdamsBashforthIntegrator(4, minStep, maxStep,
                                                                    vecAbsoluteTolerance,
                                                                    vecRelativeTolerance);
          TestProblemHandler handler = new TestProblemHandler(pb, integ);
          integ.addStepHandler(handler);
          integ.integrate(pb,
                          pb.getInitialTime(), pb.getInitialState(),
                          pb.getFinalTime(), new double[pb.getDimension()]);

    }

    @Test(expected=NumberIsTooSmallException.class)
    public void testMinStep()
            throws DimensionMismatchException, NumberIsTooSmallException,
            MaxCountExceededException, NoBracketingException {

          TestProblem1 pb = new TestProblem1();
          double minStep = 0.1 * (pb.getFinalTime() - pb.getInitialTime());
          double maxStep = pb.getFinalTime() - pb.getInitialTime();
          double[] vecAbsoluteTolerance = { 1.0e-15, 1.0e-16 };
          double[] vecRelativeTolerance = { 1.0e-15, 1.0e-16 };

          FirstOrderIntegrator integ = new AdamsMoultonIntegrator(4, minStep, maxStep,
                                                                  vecAbsoluteTolerance,
                                                                  vecRelativeTolerance);
          TestProblemHandler handler = new TestProblemHandler(pb, integ);
          integ.addStepHandler(handler);
          integ.integrate(pb,
                          pb.getInitialTime(), pb.getInitialState(),
                          pb.getFinalTime(), new double[pb.getDimension()]);

    }

[1191, 'src/test/java', 'org.apache.commons.math3.util', 'IntegerSequenceTest', 'testIncreasingRangeNegativeEnd', 69, 86]

[1191, 'src/test/java', 'org.apache.commons.math3.util', 'IntegerSequenceTest', 'testDecreasingRange', 88, 105]

    @Test
    public void testIncreasingRangeNegativeEnd() {
        final int start = -10;
        final int max = -1;
        final int step = 2;

        final List<Integer> seq = new ArrayList<Integer>();
        final IntegerSequence.Range r = IntegerSequence.range(start, max, step);
        for (Integer i : r) {
            seq.add(i);
        }

        Assert.assertEquals(5, seq.size());
        Assert.assertEquals(seq.size(), r.size());
        for (int i = 0; i < seq.size(); i++) {
            Assert.assertEquals(start + i * step, seq.get(i).intValue());
        }
    }

    @Test
    public void testDecreasingRange() {
        final int start = 10;
        final int max = -8;
        final int step = -3;

        final List<Integer> seq = new ArrayList<Integer>();
        final IntegerSequence.Range r = IntegerSequence.range(start, max, step);
        for (Integer i : r) {
            seq.add(i);
        }

        Assert.assertEquals(7, seq.size());
        Assert.assertEquals(seq.size(), r.size());
        for (int i = 0; i < seq.size(); i++) {
            Assert.assertEquals(start + i * step, seq.get(i).intValue());
        }
    }

[1200, 'src/test/java', 'org.apache.commons.math3.linear', 'BlockRealMatrixTest', 'testGetSetRowMatrixLarge', 702, 719]

[1200, 'src/test/java', 'org.apache.commons.math3.linear', 'BlockRealMatrixTest', 'testGetSetColumnMatrixLarge', 763, 781]

    @Test
    public void testGetSetRowMatrixLarge() {
        int n = 3 * BlockRealMatrix.BLOCK_SIZE;
        RealMatrix m = new BlockRealMatrix(n, n);
        RealMatrix sub = new BlockRealMatrix(1, n).scalarAdd(1);

        m.setRowMatrix(2, sub);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                if (i != 2) {
                    Assert.assertEquals(0.0, m.getEntry(i, j), 0.0);
                } else {
                    Assert.assertEquals(1.0, m.getEntry(i, j), 0.0);
                }
            }
        }
        Assert.assertEquals(sub, m.getRowMatrix(2));
    }

    @Test
    public void testGetSetColumnMatrixLarge() {
        int n = 3 * BlockRealMatrix.BLOCK_SIZE;
        RealMatrix m = new BlockRealMatrix(n, n);
        RealMatrix sub = new BlockRealMatrix(n, 1).scalarAdd(1);

        m.setColumnMatrix(2, sub);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                if (j != 2) {
                    Assert.assertEquals(0.0, m.getEntry(i, j), 0.0);
                } else {
                    Assert.assertEquals(1.0, m.getEntry(i, j), 0.0);
                }
            }
        }
        Assert.assertEquals(sub, m.getColumnMatrix(2));

    }

[1213, 'src/test/java', 'org.apache.commons.math3.random', 'RandomGeneratorAbstractTest', 'testNextInt2', 279, 293]

[1213, 'src/test/java', 'org.apache.commons.math3.random', 'RandomGeneratorAbstractTest', 'testNextLong2', 295, 309]

    @Test
    public void testNextInt2() {
        int walk = 0;
        final int N = 10000;
        for (int k = 0; k < N; ++k) {
           if (generator.nextInt() >= 0) {
               ++walk;
           } else {
               --walk;
           }
        }
        Assert.assertTrue("Walked too far astray: " + walk + "\nNote: This " +
                "test will fail randomly about 1 in 100 times.",
                FastMath.abs(walk) < FastMath.sqrt(N) * 2.576);
    }

    @Test
    public void testNextLong2() {
        int walk = 0;
        final int N = 1000;
        for (int k = 0; k < N; ++k) {
           if (generator.nextLong() >= 0) {
               ++walk;
           } else {
               --walk;
           }
        }
        Assert.assertTrue("Walked too far astray: " + walk + "\nNote: This " +
                "test will fail randomly about 1 in 100 times.",
                FastMath.abs(walk) < FastMath.sqrt(N) * 2.576);
    }

[1214, 'src/test/java', 'org.apache.commons.math3.complex', 'ComplexTest', 'testAddNaN', 108, 116]

[1214, 'src/test/java', 'org.apache.commons.math3.complex', 'ComplexTest', 'testSubtractNaN', 443, 451]

    @Test
    public void testAddNaN() {
        Complex x = new Complex(3.0, 4.0);
        Complex z = x.add(Complex.NaN);
        Assert.assertSame(Complex.NaN, z);
        z = new Complex(1, nan);
        Complex w = x.add(z);
        Assert.assertSame(Complex.NaN, w);
    }

    @Test
    public void testSubtractNaN() {
        Complex x = new Complex(3.0, 4.0);
        Complex z = x.subtract(Complex.NaN);
        Assert.assertSame(Complex.NaN, z);
        z = new Complex(1, nan);
        Complex w = x.subtract(z);
        Assert.assertSame(Complex.NaN, w);
    }

[1230, 'src/test/java', 'org.apache.commons.math3.linear', 'LUDecompositionTest', 'testSingular', 200, 209]

[1230, 'src/test/java', 'org.apache.commons.math3.linear', 'LUSolverTest', 'testSingular', 62, 71]

    @Test
    public void testSingular() {
        LUDecomposition lu =
            new LUDecomposition(MatrixUtils.createRealMatrix(testData));
        Assert.assertTrue(lu.getSolver().isNonSingular());
        lu = new LUDecomposition(MatrixUtils.createRealMatrix(singular));
        Assert.assertFalse(lu.getSolver().isNonSingular());
        lu = new LUDecomposition(MatrixUtils.createRealMatrix(bigSingular));
        Assert.assertFalse(lu.getSolver().isNonSingular());
    }

    @Test
    public void testSingular() {
        DecompositionSolver solver =
            new LUDecomposition(MatrixUtils.createRealMatrix(testData)).getSolver();
        Assert.assertTrue(solver.isNonSingular());
        solver = new LUDecomposition(MatrixUtils.createRealMatrix(singular)).getSolver();
        Assert.assertFalse(solver.isNonSingular());
        solver = new LUDecomposition(MatrixUtils.createRealMatrix(bigSingular)).getSolver();
        Assert.assertFalse(solver.isNonSingular());
    }

[1289, 'src/test/java', 'org.apache.commons.math3.geometry.euclidean.threed', 'FieldRotationDSTest', 'testArray', 1106, 1131]

[1289, 'src/test/java', 'org.apache.commons.math3.geometry.euclidean.threed', 'FieldRotationDfpTest', 'testArray', 873, 898]

    @Test
    public void testArray() throws MathIllegalArgumentException {

        FieldRotation<DerivativeStructure> r = new FieldRotation<DerivativeStructure>(createAxis(2, -3, 5),
                                                                                      createAngle(1.7),
                                                                                      RotationConvention.VECTOR_OPERATOR);

        for (double x = -0.9; x < 0.9; x += 0.2) {
            for (double y = -0.9; y < 0.9; y += 0.2) {
                for (double z = -0.9; z < 0.9; z += 0.2) {
                    FieldVector3D<DerivativeStructure> u = createVector(x, y, z);
                    FieldVector3D<DerivativeStructure> v = r.applyTo(u);
                    DerivativeStructure[] out = new DerivativeStructure[3];
                    r.applyTo(new DerivativeStructure[] { u.getX(), u.getY(), u.getZ() }, out);
                    Assert.assertEquals(v.getX().getReal(), out[0].getReal(), 1.0e-10);
                    Assert.assertEquals(v.getY().getReal(), out[1].getReal(), 1.0e-10);
                    Assert.assertEquals(v.getZ().getReal(), out[2].getReal(), 1.0e-10);
                    r.applyInverseTo(out, out);
                    Assert.assertEquals(u.getX().getReal(), out[0].getReal(), 1.0e-10);
                    Assert.assertEquals(u.getY().getReal(), out[1].getReal(), 1.0e-10);
                    Assert.assertEquals(u.getZ().getReal(), out[2].getReal(), 1.0e-10);
                }
            }
        }

    }

    @Test
    public void testArray() throws MathIllegalArgumentException {

        FieldRotation<Dfp> r = new FieldRotation<Dfp>(createAxis(2, -3, 5),
                                                      createAngle(1.7),
                                                      RotationConvention.VECTOR_OPERATOR);

        for (double x = -0.9; x < 0.9; x += 0.2) {
            for (double y = -0.9; y < 0.9; y += 0.2) {
                for (double z = -0.9; z < 0.9; z += 0.2) {
                    FieldVector3D<Dfp> u = createVector(x, y, z);
                    FieldVector3D<Dfp> v = r.applyTo(u);
                    Dfp[] out = new Dfp[3];
                    r.applyTo(new Dfp[] { u.getX(), u.getY(), u.getZ() }, out);
                    Assert.assertEquals(v.getX().getReal(), out[0].getReal(), 1.0e-10);
                    Assert.assertEquals(v.getY().getReal(), out[1].getReal(), 1.0e-10);
                    Assert.assertEquals(v.getZ().getReal(), out[2].getReal(), 1.0e-10);
                    r.applyInverseTo(out, out);
                    Assert.assertEquals(u.getX().getReal(), out[0].getReal(), 1.0e-10);
                    Assert.assertEquals(u.getY().getReal(), out[1].getReal(), 1.0e-10);
                    Assert.assertEquals(u.getZ().getReal(), out[2].getReal(), 1.0e-10);
                }
            }
        }

    }

[1305, 'src/test/java', 'org.apache.commons.math3.distribution', 'NormalDistributionTest', 'testInverseCumulativeProbabilityExtremes', 105, 111]

[1305, 'src/test/java', 'org.apache.commons.math3.distribution', 'TDistributionTest', 'testInverseCumulativeProbabilityExtremes', 96, 102]

    @Test
    public void testInverseCumulativeProbabilityExtremes() {
        setInverseCumulativeTestPoints(new double[] {0, 1});
        setInverseCumulativeTestValues(
                new double[] {Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY});
        verifyInverseCumulativeProbabilities();
    }

    @Test
    public void testInverseCumulativeProbabilityExtremes() {
        setInverseCumulativeTestPoints(new double[] {0, 1});
        setInverseCumulativeTestValues(
                new double[] {Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY});
        verifyInverseCumulativeProbabilities();
    }

[1315, 'src/test/java', 'org.apache.commons.math3.stat.inference', 'ChiSquareTestTest', 'testChiSquareZeroCount', 182, 190]

[1315, 'src/test/java', 'org.apache.commons.math3.stat.inference', 'TestUtilsTest', 'testChiSquareZeroCount', 186, 194]

    @Test
    public void testChiSquareZeroCount() {
        // Target values computed using R version 1.8.1
        long[][] counts = { {40, 0, 4}, {91, 1, 2}, {60, 2, 0}};
        Assert.assertEquals( "chi-square test statistic", 9.67444662263,
                testStatistic.chiSquare(counts), 1E-9);
        Assert.assertEquals("chi-square p-value", 0.0462835770603,
                testStatistic.chiSquareTest(counts), 1E-9);
    }

    @Test
    public void testChiSquareZeroCount() {
        // Target values computed using R version 1.8.1
        long[][] counts = { {40, 0, 4}, {91, 1, 2}, {60, 2, 0}};
        Assert.assertEquals( "chi-square test statistic", 9.67444662263,
                TestUtils.chiSquare(counts), 1E-9);
        Assert.assertEquals("chi-square p-value", 0.0462835770603,
                TestUtils.chiSquareTest(counts), 1E-9);
    }

[1317, 'src/test/java', 'org.apache.commons.math3.ode.nonstiff', 'ClassicalRungeKuttaIntegratorTest', 'testSanityChecks', 103, 131]

[1317, 'src/test/java', 'org.apache.commons.math3.ode.nonstiff', 'LutherIntegratorTest', 'testSanityChecks', 103, 131]

  @Test
  public void testSanityChecks()
      throws DimensionMismatchException, NumberIsTooSmallException,
             MaxCountExceededException, NoBracketingException {
    try  {
      TestProblem1 pb = new TestProblem1();
      new ClassicalRungeKuttaIntegrator(0.01).integrate(pb,
                                                        0.0, new double[pb.getDimension()+10],
                                                        1.0, new double[pb.getDimension()]);
        Assert.fail("an exception should have been thrown");
    } catch(DimensionMismatchException ie) {
    }
    try  {
        TestProblem1 pb = new TestProblem1();
        new ClassicalRungeKuttaIntegrator(0.01).integrate(pb,
                                                          0.0, new double[pb.getDimension()],
                                                          1.0, new double[pb.getDimension()+10]);
          Assert.fail("an exception should have been thrown");
      } catch(DimensionMismatchException ie) {
      }
    try  {
      TestProblem1 pb = new TestProblem1();
      new ClassicalRungeKuttaIntegrator(0.01).integrate(pb,
                                                        0.0, new double[pb.getDimension()],
                                                        0.0, new double[pb.getDimension()]);
        Assert.fail("an exception should have been thrown");
    } catch(NumberIsTooSmallException ie) {
    }
  }

    @Test
    public void testSanityChecks()
            throws DimensionMismatchException, NumberIsTooSmallException,
            MaxCountExceededException, NoBracketingException {
        try  {
            TestProblem1 pb = new TestProblem1();
            new LutherIntegrator(0.01).integrate(pb,
                                                 0.0, new double[pb.getDimension()+10],
                                                 1.0, new double[pb.getDimension()]);
            Assert.fail("an exception should have been thrown");
        } catch(DimensionMismatchException ie) {
        }
        try  {
            TestProblem1 pb = new TestProblem1();
            new LutherIntegrator(0.01).integrate(pb,
                                                 0.0, new double[pb.getDimension()],
                                                 1.0, new double[pb.getDimension()+10]);
            Assert.fail("an exception should have been thrown");
        } catch(DimensionMismatchException ie) {
        }
        try  {
            TestProblem1 pb = new TestProblem1();
            new LutherIntegrator(0.01).integrate(pb,
                                                 0.0, new double[pb.getDimension()],
                                                 0.0, new double[pb.getDimension()]);
            Assert.fail("an exception should have been thrown");
        } catch(NumberIsTooSmallException ie) {
        }
    }

[1336, 'src/test/java', 'org.apache.commons.math3.util', 'Decimal64Test', 'testIsInfinite', 404, 414]

[1336, 'src/test/java', 'org.apache.commons.math3.util', 'Decimal64Test', 'testIsNaN', 416, 426]

    @Test
    public void testIsInfinite() {
        Assert.assertFalse(MINUS_X.isInfinite());
        Assert.assertFalse(PLUS_X.isInfinite());
        Assert.assertFalse(MINUS_Y.isInfinite());
        Assert.assertFalse(PLUS_Y.isInfinite());
        Assert.assertFalse(Decimal64.NAN.isInfinite());

        Assert.assertTrue(Decimal64.NEGATIVE_INFINITY.isInfinite());
        Assert.assertTrue(Decimal64.POSITIVE_INFINITY.isInfinite());
    }

    @Test
    public void testIsNaN() {
        Assert.assertFalse(MINUS_X.isNaN());
        Assert.assertFalse(PLUS_X.isNaN());
        Assert.assertFalse(MINUS_Y.isNaN());
        Assert.assertFalse(PLUS_Y.isNaN());
        Assert.assertFalse(Decimal64.NEGATIVE_INFINITY.isNaN());
        Assert.assertFalse(Decimal64.POSITIVE_INFINITY.isNaN());

        Assert.assertTrue(Decimal64.NAN.isNaN());
    }

[1338, 'src/test/java', 'org.apache.commons.math3.geometry.euclidean.threed', 'Vector3DFormatAbstractTest', 'testSimpleWithDecimals', 54, 64]

[1338, 'src/test/java', 'org.apache.commons.math3.geometry.euclidean.threed', 'Vector3DFormatAbstractTest', 'testSimpleWithDecimalsTrunc', 66, 76]

    @Test
    public void testSimpleWithDecimals() {
        Vector3D c = new Vector3D(1.23, 1.43, 1.63);
        String expected =
            "{1"    + getDecimalCharacter() +
            "23; 1" + getDecimalCharacter() +
            "43; 1" + getDecimalCharacter() +
            "63}";
        String actual = vector3DFormat.format(c);
        Assert.assertEquals(expected, actual);
    }

    @Test
    public void testSimpleWithDecimalsTrunc() {
        Vector3D c = new Vector3D(1.232323232323, 1.434343434343, 1.633333333333);
        String expected =
            "{1"    + getDecimalCharacter() +
            "2323232323; 1" + getDecimalCharacter() +
            "4343434343; 1" + getDecimalCharacter() +
            "6333333333}";
        String actual = vector3DFormat.format(c);
        Assert.assertEquals(expected, actual);
    }

[1352, 'src/test/java', 'org.apache.commons.math3.linear', 'BlockFieldMatrixTest', 'testMultiply2', 344, 350]

[1352, 'src/test/java', 'org.apache.commons.math3.linear', 'FieldMatrixImplTest', 'testMultiply2', 192, 198]

    @Test
    public void testMultiply2() {
       FieldMatrix<Fraction> m3 = new BlockFieldMatrix<Fraction>(d3);
       FieldMatrix<Fraction> m4 = new BlockFieldMatrix<Fraction>(d4);
       FieldMatrix<Fraction> m5 = new BlockFieldMatrix<Fraction>(d5);
       TestUtils.assertEquals(m3.multiply(m4), m5);
   }

    @Test
    public void testMultiply2() {
       FieldMatrix<Fraction> m3 = new Array2DRowFieldMatrix<Fraction>(d3);
       FieldMatrix<Fraction> m4 = new Array2DRowFieldMatrix<Fraction>(d4);
       FieldMatrix<Fraction> m5 = new Array2DRowFieldMatrix<Fraction>(d5);
       TestUtils.assertEquals(m3.multiply(m4), m5);
   }

[1353, 'src/test/java', 'org.apache.commons.math3.fraction', 'BigFractionFormatTest', 'testFormatImproperNegative', 94, 103]

[1353, 'src/test/java', 'org.apache.commons.math3.fraction', 'FractionFormatTest', 'testFormatImproperNegative', 92, 101]

    @Test
    public void testFormatImproperNegative() {
        BigFraction c = new BigFraction(-5, 3);

        String actual = properFormat.format(c);
        Assert.assertEquals("-1 2 / 3", actual);

        actual = improperFormat.format(c);
        Assert.assertEquals("-5 / 3", actual);
    }

    @Test
    public void testFormatImproperNegative() {
        Fraction c = new Fraction(-5, 3);

        String actual = properFormat.format(c);
        Assert.assertEquals("-1 2 / 3", actual);

        actual = improperFormat.format(c);
        Assert.assertEquals("-5 / 3", actual);
    }

[1368, 'src/test/java', 'org.apache.commons.math3.fraction', 'BigFractionFormatTest', 'testParseProper', 203, 220]

[1368, 'src/test/java', 'org.apache.commons.math3.fraction', 'FractionFormatTest', 'testParseProper', 241, 258]

    @Test
    public void testParseProper() {
        String source = "1 2 / 3";

        {
            BigFraction c = properFormat.parse(source);
            Assert.assertNotNull(c);
            Assert.assertEquals(5, c.getNumeratorAsInt());
            Assert.assertEquals(3, c.getDenominatorAsInt());
        }

        try {
            improperFormat.parse(source);
            Assert.fail("invalid improper fraction.");
        } catch (MathParseException ex) {
            // success
        }
    }

    @Test
    public void testParseProper() {
        String source = "1 2 / 3";

        {
            Fraction c = properFormat.parse(source);
            Assert.assertNotNull(c);
            Assert.assertEquals(5, c.getNumerator());
            Assert.assertEquals(3, c.getDenominator());
        }

        try {
            improperFormat.parse(source);
            Assert.fail("invalid improper fraction.");
        } catch (MathParseException ex) {
            // success
        }
    }

[1379, 'src/test/java', 'org.apache.commons.math3.analysis.function', 'LogitTest', 'testDerivativesHighOrder', 97, 106]

[1379, 'src/test/java', 'org.apache.commons.math3.analysis.function', 'SigmoidTest', 'testDerivativesHighOrder', 51, 60]

    @Test
    public void testDerivativesHighOrder() {
        DerivativeStructure l = new Logit(1, 3).value(new DerivativeStructure(1, 5, 0, 1.2));
        Assert.assertEquals(-2.1972245773362193828, l.getPartialDerivative(0), 1.0e-16);
        Assert.assertEquals(5.5555555555555555555,  l.getPartialDerivative(1), 9.0e-16);
        Assert.assertEquals(-24.691358024691358025, l.getPartialDerivative(2), 2.0e-14);
        Assert.assertEquals(250.34293552812071331,  l.getPartialDerivative(3), 2.0e-13);
        Assert.assertEquals(-3749.4284407864654778, l.getPartialDerivative(4), 4.0e-12);
        Assert.assertEquals(75001.270131585632282,  l.getPartialDerivative(5), 8.0e-11);
    }

    @Test
    public void testDerivativesHighOrder() {
        DerivativeStructure s = new Sigmoid(1, 3).value(new DerivativeStructure(1, 5, 0, 1.2));
        Assert.assertEquals(2.5370495669980352859, s.getPartialDerivative(0), 5.0e-16);
        Assert.assertEquals(0.35578888129361140441, s.getPartialDerivative(1), 6.0e-17);
        Assert.assertEquals(-0.19107626464144938116,  s.getPartialDerivative(2), 6.0e-17);
        Assert.assertEquals(-0.02396830286286711696,  s.getPartialDerivative(3), 4.0e-17);
        Assert.assertEquals(0.21682059798981049049,   s.getPartialDerivative(4), 3.0e-17);
        Assert.assertEquals(-0.19186320234632658055,  s.getPartialDerivative(5), 2.0e-16);
    }

[1383, 'src/test/java', 'org.apache.commons.math3.random', 'ValueServerTest', 'testEmptyReplayFile', 126, 137]

[1383, 'src/test/java', 'org.apache.commons.math3.random', 'ValueServerTest', 'testEmptyDigestFile', 139, 150]

    @Test
    public void testEmptyReplayFile() throws Exception {
        try {
            URL url = getClass().getResource("emptyFile.txt");
            vs.setMode(ValueServer.REPLAY_MODE);
            vs.setValuesFileURL(url);
            vs.getNext();
            Assert.fail("an exception should have been thrown");
        } catch (MathIllegalStateException mise) {
            // expected behavior
        }
    }

    @Test
    public void testEmptyDigestFile() throws Exception {
        try {
            URL url = getClass().getResource("emptyFile.txt");
            vs.setMode(ValueServer.DIGEST_MODE);
            vs.setValuesFileURL(url);
            vs.computeDistribution();
            Assert.fail("an exception should have been thrown");
        } catch (ZeroException ze) {
            // expected behavior
        }
    }

[1384, 'src/test/java', 'org.apache.commons.math3.ml.neuralnet.twod', 'NeuronSquareMesh2DTest', 'testMinimalNetworkSize1', 46, 54]

[1384, 'src/test/java', 'org.apache.commons.math3.ml.neuralnet.twod', 'NeuronSquareMesh2DTest', 'testMinimalNetworkSize2', 56, 64]

    @Test(expected=NumberIsTooSmallException.class)
    public void testMinimalNetworkSize1() {
        final FeatureInitializer[] initArray = { init };

        new NeuronSquareMesh2D(1, false,
                               2, false,
                               SquareNeighbourhood.VON_NEUMANN,
                               initArray);
    }

    @Test(expected=NumberIsTooSmallException.class)
    public void testMinimalNetworkSize2() {
        final FeatureInitializer[] initArray = { init };

        new NeuronSquareMesh2D(2, false,
                               0, false,
                               SquareNeighbourhood.VON_NEUMANN,
                               initArray);
    }

[1395, 'src/test/java', 'org.apache.commons.math3.geometry.euclidean.threed', 'FieldRotationDSTest', 'testDoubleVectors', 953, 986]

[1395, 'src/test/java', 'org.apache.commons.math3.geometry.euclidean.threed', 'FieldRotationDfpTest', 'testDoubleVectors', 794, 827]

    @Test
    public void testDoubleVectors() throws MathIllegalArgumentException {

        Well1024a random = new Well1024a(0x180b41cfeeffaf67l);
        UnitSphereRandomVectorGenerator g = new UnitSphereRandomVectorGenerator(3, random);
        for (int i = 0; i < 10; ++i) {
            double[] unit = g.nextVector();
            FieldRotation<DerivativeStructure> r = new FieldRotation<DerivativeStructure>(createVector(unit[0], unit[1], unit[2]),
                                                                                          createAngle(random.nextDouble()),
                                                                                          RotationConvention.VECTOR_OPERATOR);

            for (double x = -0.9; x < 0.9; x += 0.2) {
                for (double y = -0.9; y < 0.9; y += 0.2) {
                    for (double z = -0.9; z < 0.9; z += 0.2) {
                        FieldVector3D<DerivativeStructure> uds   = createVector(x, y, z);
                        FieldVector3D<DerivativeStructure> ruds  = r.applyTo(uds);
                        FieldVector3D<DerivativeStructure> rIuds = r.applyInverseTo(uds);
                        Vector3D   u     = new Vector3D(x, y, z);
                        FieldVector3D<DerivativeStructure> ru    = r.applyTo(u);
                        FieldVector3D<DerivativeStructure> rIu   = r.applyInverseTo(u);
                        DerivativeStructure[] ruArray = new DerivativeStructure[3];
                        r.applyTo(new double[] { x, y, z}, ruArray);
                        DerivativeStructure[] rIuArray = new DerivativeStructure[3];
                        r.applyInverseTo(new double[] { x, y, z}, rIuArray);
                        checkVector(ruds, ru);
                        checkVector(ruds, new FieldVector3D<DerivativeStructure>(ruArray));
                        checkVector(rIuds, rIu);
                        checkVector(rIuds, new FieldVector3D<DerivativeStructure>(rIuArray));
                    }
                }
            }
        }

    }

    @Test
    public void testDoubleVectors() throws MathIllegalArgumentException {

        Well1024a random = new Well1024a(0x180b41cfeeffaf67l);
        UnitSphereRandomVectorGenerator g = new UnitSphereRandomVectorGenerator(3, random);
        for (int i = 0; i < 10; ++i) {
            double[] unit = g.nextVector();
            FieldRotation<Dfp> r = new FieldRotation<Dfp>(createVector(unit[0], unit[1], unit[2]),
                                                          createAngle(random.nextDouble()),
                                                          RotationConvention.VECTOR_OPERATOR);

            for (double x = -0.9; x < 0.9; x += 0.4) {
                for (double y = -0.9; y < 0.9; y += 0.4) {
                    for (double z = -0.9; z < 0.9; z += 0.4) {
                        FieldVector3D<Dfp> uds   = createVector(x, y, z);
                        FieldVector3D<Dfp> ruds  = r.applyTo(uds);
                        FieldVector3D<Dfp> rIuds = r.applyInverseTo(uds);
                        Vector3D   u     = new Vector3D(x, y, z);
                        FieldVector3D<Dfp> ru    = r.applyTo(u);
                        FieldVector3D<Dfp> rIu   = r.applyInverseTo(u);
                        Dfp[] ruArray = new Dfp[3];
                        r.applyTo(new double[] { x, y, z}, ruArray);
                        Dfp[] rIuArray = new Dfp[3];
                        r.applyInverseTo(new double[] { x, y, z}, rIuArray);
                        checkVector(ruds, ru);
                        checkVector(ruds, new FieldVector3D<Dfp>(ruArray));
                        checkVector(rIuds, rIu);
                        checkVector(rIuds, new FieldVector3D<Dfp>(rIuArray));
                    }
                }
            }
        }

    }

[1423, 'src/test/java', 'org.apache.commons.math3.linear', 'BlockFieldMatrixTest', 'testSetColumn', 1082, 1101]

[1423, 'src/test/java', 'org.apache.commons.math3.linear', 'FieldMatrixImplTest', 'testSetColumn', 821, 840]

    @Test
    public void testSetColumn() {
        FieldMatrix<Fraction> m = new BlockFieldMatrix<Fraction>(subTestData);
        Fraction[] mColumn3 = columnToArray(subColumn3);
        Assert.assertTrue(mColumn3[0] != m.getColumn(1)[0]);
        m.setColumn(1, mColumn3);
        checkArrays(mColumn3, m.getColumn(1));
        try {
            m.setColumn(-1, mColumn3);
            Assert.fail("Expecting OutOfRangeException");
        } catch (OutOfRangeException ex) {
            // expected
        }
        try {
            m.setColumn(0, new Fraction[5]);
            Assert.fail("Expecting MatrixDimensionMismatchException");
        } catch (MatrixDimensionMismatchException ex) {
            // expected
        }
    }

    @Test
    public void testSetColumn() {
        FieldMatrix<Fraction> m = new Array2DRowFieldMatrix<Fraction>(subTestData);
        Fraction[] mColumn3 = columnToArray(subColumn3);
        Assert.assertTrue(mColumn3[0] != m.getColumn(1)[0]);
        m.setColumn(1, mColumn3);
        checkArrays(mColumn3, m.getColumn(1));
        try {
            m.setColumn(-1, mColumn3);
            Assert.fail("Expecting OutOfRangeException");
        } catch (OutOfRangeException ex) {
            // expected
        }
        try {
            m.setColumn(0, new Fraction[5]);
            Assert.fail("Expecting MatrixDimensionMismatchException");
        } catch (MatrixDimensionMismatchException ex) {
            // expected
        }
    }

[1424, 'src/test/java', 'org.apache.commons.math3.complex', 'ComplexTest', 'testTanInf', 1142, 1152]

[1424, 'src/test/java', 'org.apache.commons.math3.complex', 'ComplexTest', 'testTanhInf', 1179, 1189]

    @Test
    public void testTanInf() {
        TestUtils.assertSame(Complex.valueOf(0.0, 1.0), oneInf.tan());
        TestUtils.assertSame(Complex.valueOf(0.0, -1.0), oneNegInf.tan());
        TestUtils.assertSame(Complex.NaN, infOne.tan());
        TestUtils.assertSame(Complex.NaN, negInfOne.tan());
        TestUtils.assertSame(Complex.NaN, infInf.tan());
        TestUtils.assertSame(Complex.NaN, infNegInf.tan());
        TestUtils.assertSame(Complex.NaN, negInfInf.tan());
        TestUtils.assertSame(Complex.NaN, negInfNegInf.tan());
    }

    @Test
    public void testTanhInf() {
        TestUtils.assertSame(Complex.NaN, oneInf.tanh());
        TestUtils.assertSame(Complex.NaN, oneNegInf.tanh());
        TestUtils.assertSame(Complex.valueOf(1.0, 0.0), infOne.tanh());
        TestUtils.assertSame(Complex.valueOf(-1.0, 0.0), negInfOne.tanh());
        TestUtils.assertSame(Complex.NaN, infInf.tanh());
        TestUtils.assertSame(Complex.NaN, infNegInf.tanh());
        TestUtils.assertSame(Complex.NaN, negInfInf.tanh());
        TestUtils.assertSame(Complex.NaN, negInfNegInf.tanh());
    }

[1433, 'src/test/java', 'org.apache.commons.math3.geometry.euclidean.twod.hull', 'ConvexHullGenerator2DAbstractTest', 'testOnePoint', 76, 82]

[1433, 'src/test/java', 'org.apache.commons.math3.geometry.euclidean.twod.hull', 'ConvexHullGenerator2DAbstractTest', 'testTwoPoints', 84, 90]

    @Test
    public void testOnePoint() {
        List<Vector2D> points = createRandomPoints(1);
        ConvexHull2D hull = generator.generate(points);
        Assert.assertTrue(hull.getVertices().length == 1);
        Assert.assertTrue(hull.getLineSegments().length == 0);
    }

    @Test
    public void testTwoPoints() {
        List<Vector2D> points = createRandomPoints(2);
        ConvexHull2D hull = generator.generate(points);
        Assert.assertTrue(hull.getVertices().length == 2);
        Assert.assertTrue(hull.getLineSegments().length == 1);
    }

[1446, 'src/test/java', 'org.apache.commons.math3.stat.correlation', 'KendallsCorrelationTest', 'testSimpleReversed', 93, 103]

[1446, 'src/test/java', 'org.apache.commons.math3.stat.correlation', 'KendallsCorrelationTest', 'testSimpleReversedPowerOf2', 117, 127]

    @Test
    public void testSimpleReversed() {
        final int length = 10;
        final double[] xArray = new double[length];
        final double[] yArray = new double[length];
        for (int i = 0; i < length; i++) {
            xArray[length - i - 1] = i;
            yArray[i] = i;
        }
        Assert.assertEquals(-1.0, correlation.correlation(xArray, yArray), Double.MIN_VALUE);
    }

    @Test
    public void testSimpleReversedPowerOf2() {
        final int length = 16;
        final double[] xArray = new double[length];
        final double[] yArray = new double[length];
        for (int i = 0; i < length; i++) {
            xArray[length - i - 1] = i;
            yArray[i] = i;
        }
        Assert.assertEquals(-1.0, correlation.correlation(xArray, yArray), Double.MIN_VALUE);
    }

[1453, 'src/test/java', 'org.apache.commons.math3.util', 'CombinationsTest', 'testAccessor1', 32, 37]

[1453, 'src/test/java', 'org.apache.commons.math3.util', 'CombinationsTest', 'testAccessor2', 38, 43]

    @Test
    public void testAccessor1() {
        final int n = 5;
        final int k = 3;
        Assert.assertEquals(n, new Combinations(n, k).getN());
    }

    @Test
    public void testAccessor2() {
        final int n = 5;
        final int k = 3;
        Assert.assertEquals(k, new Combinations(n, k).getK());
    }

[1472, 'src/test/java', 'org.apache.commons.math3.stat.descriptive.rank', 'PSquarePercentileTest', 'testMarkersOORLow', 171, 176]

[1472, 'src/test/java', 'org.apache.commons.math3.stat.descriptive.rank', 'PSquarePercentileTest', 'testMarkersOORHigh', 178, 183]

    @Test(expected = OutOfRangeException.class)
    public void testMarkersOORLow() {
        PSquarePercentile.newMarkers(
                Arrays.asList(new Double[] { 0.02, 1.18, 9.15, 21.91, 38.62 }),
                0.5).estimate(0);
    }

    @Test(expected = OutOfRangeException.class)
    public void testMarkersOORHigh() {
        PSquarePercentile.newMarkers(
                Arrays.asList(new Double[] { 0.02, 1.18, 9.15, 21.91, 38.62 }),
                0.5).estimate(5);
    }

[1481, 'src/test/java', 'org.apache.commons.math3.linear', 'BlockFieldMatrixTest', 'testOperateLarge', 389, 401]

[1481, 'src/test/java', 'org.apache.commons.math3.linear', 'BlockFieldMatrixTest', 'testOperatePremultiplyLarge', 403, 415]

    @Test
    public void testOperateLarge() {
        int p = (11 * BlockFieldMatrix.BLOCK_SIZE) / 10;
        int q = (11 * BlockFieldMatrix.BLOCK_SIZE) / 10;
        int r =  BlockFieldMatrix.BLOCK_SIZE / 2;
        Random random = new Random(111007463902334l);
        FieldMatrix<Fraction> m1 = createRandomMatrix(random, p, q);
        FieldMatrix<Fraction> m2 = createRandomMatrix(random, q, r);
        FieldMatrix<Fraction> m1m2 = m1.multiply(m2);
        for (int i = 0; i < r; ++i) {
            TestUtils.assertEquals(m1m2.getColumn(i), m1.operate(m2.getColumn(i)));
        }
    }

    @Test
    public void testOperatePremultiplyLarge() {
        int p = (11 * BlockFieldMatrix.BLOCK_SIZE) / 10;
        int q = (11 * BlockFieldMatrix.BLOCK_SIZE) / 10;
        int r =  BlockFieldMatrix.BLOCK_SIZE / 2;
        Random random = new Random(111007463902334l);
        FieldMatrix<Fraction> m1 = createRandomMatrix(random, p, q);
        FieldMatrix<Fraction> m2 = createRandomMatrix(random, q, r);
        FieldMatrix<Fraction> m1m2 = m1.multiply(m2);
        for (int i = 0; i < p; ++i) {
            TestUtils.assertEquals(m1m2.getRow(i), m2.preMultiply(m1.getRow(i)));
        }
    }

[1562, 'src/test/java', 'org.apache.commons.math3.optim.nonlinear.scalar.noderiv', 'SimplexOptimizerMultiDirectionalTest', 'testBoundsUnsupported', 34, 46]

[1562, 'src/test/java', 'org.apache.commons.math3.optim.nonlinear.scalar.noderiv', 'SimplexOptimizerNelderMeadTest', 'testBoundsUnsupported', 39, 51]

    @Test(expected=MathUnsupportedOperationException.class)
    public void testBoundsUnsupported() {
        SimplexOptimizer optimizer = new SimplexOptimizer(1e-10, 1e-30);
        final FourExtrema fourExtrema = new FourExtrema();

        optimizer.optimize(new MaxEval(100),
                           new ObjectiveFunction(fourExtrema),
                           GoalType.MINIMIZE,
                           new InitialGuess(new double[] { -3, 0 }),
                           new NelderMeadSimplex(new double[] { 0.2, 0.2 }),
                           new SimpleBounds(new double[] { -5, -1 },
                                            new double[] { 5, 1 }));
    }

    @Test(expected=MathUnsupportedOperationException.class)
    public void testBoundsUnsupported() {
        SimplexOptimizer optimizer = new SimplexOptimizer(1e-10, 1e-30);
        final FourExtrema fourExtrema = new FourExtrema();

        optimizer.optimize(new MaxEval(100),
                           new ObjectiveFunction(fourExtrema),
                           GoalType.MINIMIZE,
                           new InitialGuess(new double[] { -3, 0 }),
                           new NelderMeadSimplex(new double[] { 0.2, 0.2 }),
                           new SimpleBounds(new double[] { -5, -1 },
                                            new double[] { 5, 1 }));
    }

[1564, 'src/test/java', 'org.apache.commons.math3.stat.inference', 'ChiSquareTestTest', 'testChiSquareLargeTestStatistic', 163, 179]

[1564, 'src/test/java', 'org.apache.commons.math3.stat.inference', 'TestUtilsTest', 'testChiSquareLargeTestStatistic', 167, 183]

    @Test
    public void testChiSquareLargeTestStatistic() {
        double[] exp = new double[] {
            3389119.5, 649136.6, 285745.4, 25357364.76, 11291189.78, 543628.0,
            232921.0, 437665.75
        };

        long[] obs = new long[] {
            2372383, 584222, 257170, 17750155, 7903832, 489265, 209628, 393899
        };
        org.apache.commons.math3.stat.inference.ChiSquareTest csti =
            new org.apache.commons.math3.stat.inference.ChiSquareTest();
        double cst = csti.chiSquareTest(exp, obs);
        Assert.assertEquals("chi-square p-value", 0.0, cst, 1E-3);
        Assert.assertEquals( "chi-square test statistic",
                114875.90421929007, testStatistic.chiSquare(exp, obs), 1E-9);
    }

    @Test
    public void testChiSquareLargeTestStatistic() {
        double[] exp = new double[] {
                3389119.5, 649136.6, 285745.4, 25357364.76, 11291189.78, 543628.0,
                232921.0, 437665.75
        };

        long[] obs = new long[] {
                2372383, 584222, 257170, 17750155, 7903832, 489265, 209628, 393899
        };
        org.apache.commons.math3.stat.inference.ChiSquareTest csti =
            new org.apache.commons.math3.stat.inference.ChiSquareTest();
        double cst = csti.chiSquareTest(exp, obs);
        Assert.assertEquals("chi-square p-value", 0.0, cst, 1E-3);
        Assert.assertEquals( "chi-square test statistic",
                114875.90421929007, TestUtils.chiSquare(exp, obs), 1E-9);
    }

[1573, 'src/test/java', 'org.apache.commons.math3.linear', 'RealVectorAbstractTest', 'testSetSubVectorInvalidIndex2', 430, 433]

[1573, 'src/test/java', 'org.apache.commons.math3.linear', 'RealVectorAbstractTest', 'testSetSubVectorInvalidIndex3', 435, 438]

    @Test(expected = OutOfRangeException.class)
    public void testSetSubVectorInvalidIndex2() {
        create(new double[10]).setSubVector(10, create(new double[2]));
    }

    @Test(expected = OutOfRangeException.class)
    public void testSetSubVectorInvalidIndex3() {
        create(new double[10]).setSubVector(9, create(new double[2]));
    }

[1576, 'src/test/java', 'org.apache.commons.math3.random', 'StableRandomGeneratorTest', 'testAlphaRangeBelowZero', 88, 97]

[1576, 'src/test/java', 'org.apache.commons.math3.random', 'StableRandomGeneratorTest', 'testBetaRangeBelowMinusOne', 110, 119]

    @Test
    public void testAlphaRangeBelowZero() {
        try {
            new StableRandomGenerator(rg,
                    -1.0, 0.0);
            Assert.fail("Expected OutOfRangeException");
        } catch (OutOfRangeException e) {
            Assert.assertEquals(-1.0, e.getArgument());
        }
    }

    @Test
    public void testBetaRangeBelowMinusOne() {
        try {
            new StableRandomGenerator(rg,
                    1.0, -2.0);
            Assert.fail("Expected OutOfRangeException");
        } catch (OutOfRangeException e) {
            Assert.assertEquals(-2.0, e.getArgument());
        }
    }

[1588, 'src/test/java', 'org.apache.commons.math3.analysis.differentiation', 'SparseGradientTest', 'testToDegreesDefinition', 910, 920]

[1588, 'src/test/java', 'org.apache.commons.math3.analysis.differentiation', 'SparseGradientTest', 'testToRadiansDefinition', 922, 932]

    @Test
    public void testToDegreesDefinition() {
        double epsilon = 3.0e-16;
        for (int maxOrder = 0; maxOrder < 6; ++maxOrder) {
            for (double x = 0.1; x < 1.2; x += 0.001) {
                SparseGradient sgX = SparseGradient.createVariable(0, x);
                Assert.assertEquals(FastMath.toDegrees(x), sgX.toDegrees().getValue(), epsilon);
                Assert.assertEquals(180 / FastMath.PI, sgX.toDegrees().getDerivative(0), epsilon);
            }
        }
    }

    @Test
    public void testToRadiansDefinition() {
        double epsilon = 3.0e-16;
        for (int maxOrder = 0; maxOrder < 6; ++maxOrder) {
            for (double x = 0.1; x < 1.2; x += 0.001) {
                SparseGradient sgX = SparseGradient.createVariable(0, x);
                Assert.assertEquals(FastMath.toRadians(x), sgX.toRadians().getValue(), epsilon);
                Assert.assertEquals(FastMath.PI / 180, sgX.toRadians().getDerivative(0), epsilon);
            }
        }
    }

[1593, 'src/test/java', 'org.apache.commons.math3.linear', 'BlockFieldMatrixTest', 'testEqualsAndHashCode', 1139, 1151]

[1593, 'src/test/java', 'org.apache.commons.math3.linear', 'FieldMatrixImplTest', 'testEqualsAndHashCode', 857, 869]

    @Test
    public void testEqualsAndHashCode() {
        BlockFieldMatrix<Fraction> m = new BlockFieldMatrix<Fraction>(testData);
        BlockFieldMatrix<Fraction> m1 = (BlockFieldMatrix<Fraction>) m.copy();
        BlockFieldMatrix<Fraction> mt = (BlockFieldMatrix<Fraction>) m.transpose();
        Assert.assertTrue(m.hashCode() != mt.hashCode());
        Assert.assertEquals(m.hashCode(), m1.hashCode());
        Assert.assertEquals(m, m);
        Assert.assertEquals(m, m1);
        Assert.assertFalse(m.equals(null));
        Assert.assertFalse(m.equals(mt));
        Assert.assertFalse(m.equals(new BlockFieldMatrix<Fraction>(bigSingular)));
    }

    @Test
    public void testEqualsAndHashCode() {
        Array2DRowFieldMatrix<Fraction> m = new Array2DRowFieldMatrix<Fraction>(testData);
        Array2DRowFieldMatrix<Fraction> m1 = (Array2DRowFieldMatrix<Fraction>) m.copy();
        Array2DRowFieldMatrix<Fraction> mt = (Array2DRowFieldMatrix<Fraction>) m.transpose();
        Assert.assertTrue(m.hashCode() != mt.hashCode());
        Assert.assertEquals(m.hashCode(), m1.hashCode());
        Assert.assertEquals(m, m);
        Assert.assertEquals(m, m1);
        Assert.assertFalse(m.equals(null));
        Assert.assertFalse(m.equals(mt));
        Assert.assertFalse(m.equals(new Array2DRowFieldMatrix<Fraction>(bigSingular)));
    }

[1595, 'src/test/java', 'org.apache.commons.math3.stat.inference', 'ChiSquareTestTest', 'testChiSquare', 37, 102]

[1595, 'src/test/java', 'org.apache.commons.math3.stat.inference', 'TestUtilsTest', 'testChiSquare', 41, 106]

    @Test
    public void testChiSquare() {

        // Target values computed using R version 1.8.1
        // Some assembly required ;-)
        //      Use sum((obs - exp)^2/exp) for the chi-square statistic and
        //      1 - pchisq(sum((obs - exp)^2/exp), length(obs) - 1) for the p-value

        long[] observed = {10, 9, 11};
        double[] expected = {10, 10, 10};
        Assert.assertEquals("chi-square statistic", 0.2,  testStatistic.chiSquare(expected, observed), 10E-12);
        Assert.assertEquals("chi-square p-value", 0.904837418036, testStatistic.chiSquareTest(expected, observed), 1E-10);

        long[] observed1 = { 500, 623, 72, 70, 31 };
        double[] expected1 = { 485, 541, 82, 61, 37 };
        Assert.assertEquals( "chi-square test statistic", 9.023307936427388, testStatistic.chiSquare(expected1, observed1), 1E-10);
        Assert.assertEquals("chi-square p-value", 0.06051952647453607, testStatistic.chiSquareTest(expected1, observed1), 1E-9);
        Assert.assertTrue("chi-square test reject", testStatistic.chiSquareTest(expected1, observed1, 0.08));
        Assert.assertTrue("chi-square test accept", !testStatistic.chiSquareTest(expected1, observed1, 0.05));

        try {
            testStatistic.chiSquareTest(expected1, observed1, 95);
            Assert.fail("alpha out of range, OutOfRangeException expected");
        } catch (OutOfRangeException ex) {
            // expected
        }

        long[] tooShortObs = { 0 };
        double[] tooShortEx = { 1 };
        try {
            testStatistic.chiSquare(tooShortEx, tooShortObs);
            Assert.fail("arguments too short, DimensionMismatchException expected");
        } catch (DimensionMismatchException ex) {
            // expected
        }

        // unmatched arrays
        long[] unMatchedObs = { 0, 1, 2, 3 };
        double[] unMatchedEx = { 1, 1, 2 };
        try {
            testStatistic.chiSquare(unMatchedEx, unMatchedObs);
            Assert.fail("arrays have different lengths, DimensionMismatchException expected");
        } catch (DimensionMismatchException ex) {
            // expected
        }

        // 0 expected count
        expected[0] = 0;
        try {
            testStatistic.chiSquareTest(expected, observed, .01);
            Assert.fail("bad expected count, NotStrictlyPositiveException expected");
        } catch (NotStrictlyPositiveException ex) {
            // expected
        }

        // negative observed count
        expected[0] = 1;
        observed[0] = -1;
        try {
            testStatistic.chiSquareTest(expected, observed, .01);
            Assert.fail("bad expected count, NotPositiveException expected");
        } catch (NotPositiveException ex) {
            // expected
        }

    }

    @Test
    public void testChiSquare() {

        // Target values computed using R version 1.8.1
        // Some assembly required ;-)
        //      Use sum((obs - exp)^2/exp) for the chi-square statistic and
        //      1 - pchisq(sum((obs - exp)^2/exp), length(obs) - 1) for the p-value

        long[] observed = {10, 9, 11};
        double[] expected = {10, 10, 10};
        Assert.assertEquals("chi-square statistic", 0.2,  TestUtils.chiSquare(expected, observed), 10E-12);
        Assert.assertEquals("chi-square p-value", 0.904837418036, TestUtils.chiSquareTest(expected, observed), 1E-10);

        long[] observed1 = { 500, 623, 72, 70, 31 };
        double[] expected1 = { 485, 541, 82, 61, 37 };
        Assert.assertEquals( "chi-square test statistic", 9.023307936427388, TestUtils.chiSquare(expected1, observed1), 1E-10);
        Assert.assertEquals("chi-square p-value", 0.06051952647453607, TestUtils.chiSquareTest(expected1, observed1), 1E-9);
        Assert.assertTrue("chi-square test reject", TestUtils.chiSquareTest(expected1, observed1, 0.07));
        Assert.assertTrue("chi-square test accept", !TestUtils.chiSquareTest(expected1, observed1, 0.05));

        try {
            TestUtils.chiSquareTest(expected1, observed1, 95);
            Assert.fail("alpha out of range, OutOfRangeException expected");
        } catch (OutOfRangeException ex) {
            // expected
        }

        long[] tooShortObs = { 0 };
        double[] tooShortEx = { 1 };
        try {
            TestUtils.chiSquare(tooShortEx, tooShortObs);
            Assert.fail("arguments too short, DimensionMismatchException expected");
        } catch (DimensionMismatchException ex) {
            // expected
        }

        // unmatched arrays
        long[] unMatchedObs = { 0, 1, 2, 3 };
        double[] unMatchedEx = { 1, 1, 2 };
        try {
            TestUtils.chiSquare(unMatchedEx, unMatchedObs);
            Assert.fail("arrays have different lengths, DimensionMismatchException expected");
        } catch (DimensionMismatchException ex) {
            // expected
        }

        // 0 expected count
        expected[0] = 0;
        try {
            TestUtils.chiSquareTest(expected, observed, .01);
            Assert.fail("bad expected count, NotStrictlyPositiveException expected");
        } catch (NotStrictlyPositiveException ex) {
            // expected
        }

        // negative observed count
        expected[0] = 1;
        observed[0] = -1;
        try {
            TestUtils.chiSquareTest(expected, observed, .01);
            Assert.fail("bad expected count, NotPositiveException expected");
        } catch (NotPositiveException ex) {
            // expected
        }

    }

[1611, 'src/test/java', 'org.apache.commons.math3.distribution', 'FDistributionTest', 'testPreconditions', 90, 104]

[1611, 'src/test/java', 'org.apache.commons.math3.distribution', 'GammaDistributionTest', 'testPreconditions', 85, 99]

    @Test
    public void testPreconditions() {
        try {
            new FDistribution(0, 1);
            Assert.fail("Expecting NotStrictlyPositiveException for df = 0");
        } catch (NotStrictlyPositiveException ex) {
            // Expected.
        }
        try {
            new FDistribution(1, 0);
            Assert.fail("Expecting NotStrictlyPositiveException for df = 0");
        } catch (NotStrictlyPositiveException ex) {
            // Expected.
        }
    }

    @Test
    public void testPreconditions() {
        try {
            new GammaDistribution(0, 1);
            Assert.fail("Expecting NotStrictlyPositiveException for alpha = 0");
        } catch (NotStrictlyPositiveException ex) {
            // Expected.
        }
        try {
            new GammaDistribution(1, 0);
            Assert.fail("Expecting NotStrictlyPositiveException for alpha = 0");
        } catch (NotStrictlyPositiveException ex) {
            // Expected.
        }
    }

[1637, 'src/test/java', 'org.apache.commons.math3.linear', 'RealMatrixFormatAbstractTest', 'testParseNoComponents', 316, 324]

[1637, 'src/test/java', 'org.apache.commons.math3.linear', 'RealVectorFormatAbstractTest', 'testParseNoComponents', 308, 316]

    @Test
    public void testParseNoComponents() {
        try {
            realMatrixFormat.parse("{{ }}");
            Assert.fail("Expecting MathParseException");
        } catch (MathParseException pe) {
            // expected behavior
        }
    }

    @Test
    public void testParseNoComponents() {
        try {
            realVectorFormat.parse("{ }");
            Assert.fail("Expecting MathParseException");
        } catch (MathParseException pe) {
            // expected behavior
        }
    }

[1659, 'src/test/java', 'org.apache.commons.math3.analysis.integration.gauss', 'GaussIntegratorTest', 'testGetWeights', 30, 43]

[1659, 'src/test/java', 'org.apache.commons.math3.analysis.integration.gauss', 'GaussIntegratorTest', 'testGetPoints', 45, 58]

    @Test
    public void testGetWeights() {
        final double[] points = { 0, 1.2, 3.4 };
        final double[] weights = { 9.8, 7.6, 5.4 };

        final GaussIntegrator integrator
            = new GaussIntegrator(new Pair<double[], double[]>(points, weights));

        Assert.assertEquals(weights.length, integrator.getNumberOfPoints());

        for (int i = 0; i < integrator.getNumberOfPoints(); i++) {
            Assert.assertEquals(weights[i], integrator.getWeight(i), 0d);
        }
    }

    @Test
    public void testGetPoints() {
        final double[] points = { 0, 1.2, 3.4 };
        final double[] weights = { 9.8, 7.6, 5.4 };

        final GaussIntegrator integrator
            = new GaussIntegrator(new Pair<double[], double[]>(points, weights));

        Assert.assertEquals(points.length, integrator.getNumberOfPoints());

        for (int i = 0; i < integrator.getNumberOfPoints(); i++) {
            Assert.assertEquals(points[i], integrator.getPoint(i), 0d);
        }
    }

[1682, 'src/test/java', 'org.apache.commons.math3.stat.regression', 'SimpleRegressionTest', 'testRemoveSingle', 621, 637]

[1682, 'src/test/java', 'org.apache.commons.math3.stat.regression', 'SimpleRegressionTest', 'testRemoveMultiple', 640, 656]

    @Test
    public void testRemoveSingle() {
        // Create regression with inference data then remove to test
        SimpleRegression regression = new SimpleRegression();
        regression.addData(infData);
        regression.removeData(removeSingle);
        regression.addData(removeSingle);
        // Use the inference assertions to make sure that everything worked
        Assert.assertEquals("slope std err", 0.011448491,
                regression.getSlopeStdErr(), 1E-10);
        Assert.assertEquals("std err intercept", 0.286036932,
                regression.getInterceptStdErr(),1E-8);
        Assert.assertEquals("significance", 4.596e-07,
                regression.getSignificance(),1E-8);
        Assert.assertEquals("slope conf interval half-width", 0.0270713794287,
                regression.getSlopeConfidenceInterval(),1E-8);
     }

    @Test
    public void testRemoveMultiple() {
        // Create regression with inference data then remove to test
        SimpleRegression regression = new SimpleRegression();
        regression.addData(infData);
        regression.removeData(removeMultiple);
        regression.addData(removeMultiple);
        // Use the inference assertions to make sure that everything worked
        Assert.assertEquals("slope std err", 0.011448491,
                regression.getSlopeStdErr(), 1E-10);
        Assert.assertEquals("std err intercept", 0.286036932,
                regression.getInterceptStdErr(),1E-8);
        Assert.assertEquals("significance", 4.596e-07,
                regression.getSignificance(),1E-8);
        Assert.assertEquals("slope conf interval half-width", 0.0270713794287,
                regression.getSlopeConfidenceInterval(),1E-8);
     }

[1686, 'src/test/java', 'org.apache.commons.math3.linear', 'Array2DRowRealMatrixTest', 'testAdd', 125, 138]

[1686, 'src/test/java', 'org.apache.commons.math3.linear', 'BlockRealMatrixTest', 'testAdd', 126, 139]

    @Test
    public void testAdd() {
        Array2DRowRealMatrix m = new Array2DRowRealMatrix(testData);
        Array2DRowRealMatrix mInv = new Array2DRowRealMatrix(testDataInv);
        RealMatrix mPlusMInv = m.add(mInv);
        double[][] sumEntries = mPlusMInv.getData();
        for (int row = 0; row < m.getRowDimension(); row++) {
            for (int col = 0; col < m.getColumnDimension(); col++) {
                Assert.assertEquals("sum entry entry",
                    testDataPlusInv[row][col],sumEntries[row][col],
                        entryTolerance);
            }
        }
    }

    @Test
    public void testAdd() {
        BlockRealMatrix m = new BlockRealMatrix(testData);
        BlockRealMatrix mInv = new BlockRealMatrix(testDataInv);
        RealMatrix mPlusMInv = m.add(mInv);
        double[][] sumEntries = mPlusMInv.getData();
        for (int row = 0; row < m.getRowDimension(); row++) {
            for (int col = 0; col < m.getColumnDimension(); col++) {
                Assert.assertEquals("sum entry entry",
                    testDataPlusInv[row][col],sumEntries[row][col],
                        entryTolerance);
            }
        }
    }

[1695, 'src/test/java', 'org.apache.commons.math3.geometry.euclidean.twod', 'PolygonsSetTest', 'testEmpty', 119, 133]

[1695, 'src/test/java', 'org.apache.commons.math3.geometry.euclidean.twod', 'PolygonsSetTest', 'testFull', 135, 149]

    @Test
    public void testEmpty() {
        PolygonsSet empty = (PolygonsSet) new RegionFactory<Euclidean2D>().getComplement(new PolygonsSet(1.0e-10));
        Assert.assertTrue(empty.isEmpty());
        Assert.assertEquals(0, empty.getVertices().length);
        Assert.assertEquals(0.0, empty.getBoundarySize(), 1.0e-10);
        Assert.assertEquals(0.0, empty.getSize(), 1.0e-10);
        for (double y = -1; y < 1; y += 0.1) {
            for (double x = -1; x < 1; x += 0.1) {
                Assert.assertEquals(Double.POSITIVE_INFINITY,
                                    empty.projectToBoundary(new Vector2D(x, y)).getOffset(),
                                    1.0e-10);
            }
        }
    }

    @Test
    public void testFull() {
        PolygonsSet empty = new PolygonsSet(1.0e-10);
        Assert.assertFalse(empty.isEmpty());
        Assert.assertEquals(0, empty.getVertices().length);
        Assert.assertEquals(0.0, empty.getBoundarySize(), 1.0e-10);
        Assert.assertEquals(Double.POSITIVE_INFINITY, empty.getSize(), 1.0e-10);
        for (double y = -1; y < 1; y += 0.1) {
            for (double x = -1; x < 1; x += 0.1) {
                Assert.assertEquals(Double.NEGATIVE_INFINITY,
                                    empty.projectToBoundary(new Vector2D(x, y)).getOffset(),
                                    1.0e-10);
            }
        }
    }

[1719, 'src/test/java', 'org.apache.commons.math3.linear', 'RealVectorAbstractTest', 'testSetEntryInvalidIndex1', 252, 255]

[1719, 'src/test/java', 'org.apache.commons.math3.linear', 'RealVectorAbstractTest', 'testAddToEntryInvalidIndex1', 299, 302]

    @Test(expected=OutOfRangeException.class)
    public void testSetEntryInvalidIndex1() {
        create(new double[4]).setEntry(-1, getPreferredEntryValue());
    }

    @Test(expected=OutOfRangeException.class)
    public void testAddToEntryInvalidIndex1() {
        create(new double[3]).addToEntry(-1, getPreferredEntryValue());
    }

[1721, 'src/test/java', 'org.apache.commons.math3.linear', 'EigenSolverTest', 'testNonInvertibleMath1045', 96, 101]

[1721, 'src/test/java', 'org.apache.commons.math3.linear', 'QRDecompositionTest', 'testNonInvertible', 244, 249]

    @Test(expected=SingularMatrixException.class)
    public void testNonInvertibleMath1045() {
        EigenDecomposition eigen =
            new EigenDecomposition(MatrixUtils.createRealMatrix(bigSingular));
        eigen.getSolver().getInverse();
    }

    @Test(expected=SingularMatrixException.class)
    public void testNonInvertible() {
        QRDecomposition qr =
            new QRDecomposition(MatrixUtils.createRealMatrix(testData3x3Singular));
        qr.getSolver().getInverse();
    }

[1732, 'src/test/java', 'org.apache.commons.math3.linear', 'MatrixUtilsTest', 'testInverseSingular', 451, 455]

[1732, 'src/test/java', 'org.apache.commons.math3.linear', 'MatrixUtilsTest', 'testInverseNonSquare', 457, 461]

    @Test(expected=SingularMatrixException.class)
    public void testInverseSingular() {
        RealMatrix m = MatrixUtils.createRealMatrix(testData3x3Singular);
        MatrixUtils.inverse(m);
    }

    @Test(expected=NonSquareMatrixException.class)
    public void testInverseNonSquare() {
        RealMatrix m = MatrixUtils.createRealMatrix(testData3x4);
        MatrixUtils.inverse(m);
    }

[1762, 'src/test/java', 'org.apache.commons.math3.special', 'GammaTest', 'testDigammaNonRealArgs', 127, 132]

[1762, 'src/test/java', 'org.apache.commons.math3.special', 'GammaTest', 'testTrigammaNonRealArgs', 160, 165]

    @Test
    public void testDigammaNonRealArgs() {
        Assert.assertTrue(Double.isNaN(Gamma.digamma(Double.NaN)));
        Assert.assertTrue(Double.isInfinite(Gamma.digamma(Double.POSITIVE_INFINITY)));
        Assert.assertTrue(Double.isInfinite(Gamma.digamma(Double.NEGATIVE_INFINITY)));
    }

    @Test
    public void testTrigammaNonRealArgs() {
        Assert.assertTrue(Double.isNaN(Gamma.trigamma(Double.NaN)));
        Assert.assertTrue(Double.isInfinite(Gamma.trigamma(Double.POSITIVE_INFINITY)));
        Assert.assertTrue(Double.isInfinite(Gamma.trigamma(Double.NEGATIVE_INFINITY)));
    }

[1763, 'src/test/java', 'org.apache.commons.math3.stat.correlation', 'KendallsCorrelationTest', 'testAllTiesInBoth', 173, 179]

[1763, 'src/test/java', 'org.apache.commons.math3.stat.correlation', 'KendallsCorrelationTest', 'testSingleElement', 203, 209]

    @Test
    public void testAllTiesInBoth() {
        final int length = 10;
        final double[] xArray = new double[length];
        final double[] yArray = new double[length];
        Assert.assertEquals(Double.NaN, correlation.correlation(xArray, yArray), 0);
    }

    @Test
    public void testSingleElement() {
        final int length = 1;
        final double[] xArray = new double[length];
        final double[] yArray = new double[length];
        Assert.assertEquals(Double.NaN, correlation.correlation(xArray, yArray), 0);
    }

[1781, 'src/test/java', 'org.apache.commons.math3.random', 'StableRandomGeneratorTest', 'testAlphaRangeAboveTwo', 99, 108]

[1781, 'src/test/java', 'org.apache.commons.math3.random', 'StableRandomGeneratorTest', 'testBetaRangeAboveOne', 121, 130]

    @Test
    public void testAlphaRangeAboveTwo() {
        try {
            new StableRandomGenerator(rg,
                    3.0, 0.0);
            Assert.fail("Expected OutOfRangeException");
        } catch (OutOfRangeException e) {
            Assert.assertEquals(3.0, e.getArgument());
        }
    }

    @Test
    public void testBetaRangeAboveOne() {
        try {
            new StableRandomGenerator(rg,
                    1.0, 2.0);
            Assert.fail("Expected OutOfRangeException");
        } catch (OutOfRangeException e) {
            Assert.assertEquals(2.0, e.getArgument());
        }
    }

[1782, 'src/test/java', 'org.apache.commons.math3.random', 'RandomDataGeneratorTest', 'testNextIntNegativeToPositiveRange', 108, 114]

[1782, 'src/test/java', 'org.apache.commons.math3.random', 'RandomDataGeneratorTest', 'testNextLongNegativeToPositiveRange', 183, 189]

    @Test
    public void testNextIntNegativeToPositiveRange() {
        for (int i = 0; i < 5; i++) {
            checkNextIntUniform(-3, 5);
            checkNextIntUniform(-3, 6);
        }
    }

    @Test
    public void testNextLongNegativeToPositiveRange() {
        for (int i = 0; i < 5; i++) {
            checkNextLongUniform(-3, 5);
            checkNextLongUniform(-3, 6);
        }
    }

[1800, 'src/test/java', 'org.apache.commons.math3.stat.inference', 'TTestTest', 'testPaired', 284, 296]

[1800, 'src/test/java', 'org.apache.commons.math3.stat.inference', 'TestUtilsTest', 'testPaired', 437, 449]

    @Test
    public void testPaired() {
        double[] sample1 = {1d, 3d, 5d, 7d};
        double[] sample2 = {0d, 6d, 11d, 2d};
        double[] sample3 = {5d, 7d, 8d, 10d};

        // Target values computed using R, version 1.8.1 (linux version)
        Assert.assertEquals(-0.3133, testStatistic.pairedT(sample1, sample2), 1E-4);
        Assert.assertEquals(0.774544295819, testStatistic.pairedTTest(sample1, sample2), 1E-10);
        Assert.assertEquals(0.001208, testStatistic.pairedTTest(sample1, sample3), 1E-6);
        Assert.assertFalse(testStatistic.pairedTTest(sample1, sample3, .001));
        Assert.assertTrue(testStatistic.pairedTTest(sample1, sample3, .002));
    }

    @Test
    public void testPaired() {
        double[] sample1 = {1d, 3d, 5d, 7d};
        double[] sample2 = {0d, 6d, 11d, 2d};
        double[] sample3 = {5d, 7d, 8d, 10d};

        // Target values computed using R, version 1.8.1 (linux version)
        Assert.assertEquals(-0.3133, TestUtils.pairedT(sample1, sample2), 1E-4);
        Assert.assertEquals(0.774544295819, TestUtils.pairedTTest(sample1, sample2), 1E-10);
        Assert.assertEquals(0.001208, TestUtils.pairedTTest(sample1, sample3), 1E-6);
        Assert.assertFalse(TestUtils.pairedTTest(sample1, sample3, .001));
        Assert.assertTrue(TestUtils.pairedTTest(sample1, sample3, .002));
    }

[1804, 'src/test/java', 'org.apache.commons.math3.linear', 'Array2DRowRealMatrixTest', 'testSetRowVector', 723, 742]

[1804, 'src/test/java', 'org.apache.commons.math3.linear', 'BlockRealMatrixTest', 'testSetRowVector', 804, 823]

    @Test
    public void testSetRowVector() {
        RealMatrix m = new Array2DRowRealMatrix(subTestData);
        RealVector mRow3 = new ArrayRealVector(subRow3[0]);
        Assert.assertNotSame(mRow3, m.getRowMatrix(0));
        m.setRowVector(0, mRow3);
        Assert.assertEquals(mRow3, m.getRowVector(0));
        try {
            m.setRowVector(-1, mRow3);
            Assert.fail("Expecting OutOfRangeException");
        } catch (OutOfRangeException ex) {
            // expected
        }
        try {
            m.setRowVector(0, new ArrayRealVector(5));
            Assert.fail("Expecting MatrixDimensionMismatchException");
        } catch (MatrixDimensionMismatchException ex) {
            // expected
        }
    }

    @Test
    public void testSetRowVector() {
        RealMatrix m = new BlockRealMatrix(subTestData);
        RealVector mRow3 = new ArrayRealVector(subRow3[0]);
        Assert.assertNotSame(mRow3, m.getRowMatrix(0));
        m.setRowVector(0, mRow3);
        Assert.assertEquals(mRow3, m.getRowVector(0));
        try {
            m.setRowVector(-1, mRow3);
            Assert.fail("Expecting OutOfRangeException");
        } catch (OutOfRangeException ex) {
            // expected
        }
        try {
            m.setRowVector(0, new ArrayRealVector(5));
            Assert.fail("Expecting MatrixDimensionMismatchException");
        } catch (MatrixDimensionMismatchException ex) {
            // expected
        }
    }

[1810, 'src/test/java', 'org.apache.commons.math3.linear', 'ConjugateGradientTest', 'testUnpreconditionedNormOfResidual', 510, 554]

[1810, 'src/test/java', 'org.apache.commons.math3.linear', 'SymmLQTest', 'testUnpreconditionedNormOfResidual', 594, 638]

    @Test
    public void testUnpreconditionedNormOfResidual() {
        final int n = 5;
        final int maxIterations = 100;
        final RealLinearOperator a = new HilbertMatrix(n);
        final IterativeLinearSolver solver;
        final IterationListener listener = new IterationListener() {

            private void doTestNormOfResidual(final IterationEvent e) {
                final IterativeLinearSolverEvent evt;
                evt = (IterativeLinearSolverEvent) e;
                final RealVector x = evt.getSolution();
                final RealVector b = evt.getRightHandSideVector();
                final RealVector r = b.subtract(a.operate(x));
                final double rnorm = r.getNorm();
                Assert.assertEquals("iteration performed (residual)",
                    rnorm, evt.getNormOfResidual(),
                    FastMath.max(1E-5 * rnorm, 1E-10));
            }

            public void initializationPerformed(final IterationEvent e) {
                doTestNormOfResidual(e);
            }

            public void iterationPerformed(final IterationEvent e) {
                doTestNormOfResidual(e);
            }

            public void iterationStarted(final IterationEvent e) {
                doTestNormOfResidual(e);
            }

            public void terminationPerformed(final IterationEvent e) {
                doTestNormOfResidual(e);
            }
        };
        solver = new ConjugateGradient(maxIterations, 1E-10, true);
        solver.getIterationManager().addIterationListener(listener);
        final RealVector b = new ArrayRealVector(n);
        for (int j = 0; j < n; j++) {
            b.set(0.);
            b.setEntry(j, 1.);
            solver.solve(a, b);
        }
    }

    @Test
    public void testUnpreconditionedNormOfResidual() {
        final int n = 5;
        final int maxIterations = 100;
        final RealLinearOperator a = new HilbertMatrix(n);
        final IterativeLinearSolver solver;
        final IterationListener listener = new IterationListener() {

            private void doTestNormOfResidual(final IterationEvent e) {
                final IterativeLinearSolverEvent evt;
                evt = (IterativeLinearSolverEvent) e;
                final RealVector x = evt.getSolution();
                final RealVector b = evt.getRightHandSideVector();
                final RealVector r = b.subtract(a.operate(x));
                final double rnorm = r.getNorm();
                Assert.assertEquals("iteration performed (residual)",
                    rnorm, evt.getNormOfResidual(),
                    FastMath.max(1E-5 * rnorm, 1E-10));
            }

            public void initializationPerformed(final IterationEvent e) {
                doTestNormOfResidual(e);
            }

            public void iterationPerformed(final IterationEvent e) {
                doTestNormOfResidual(e);
            }

            public void iterationStarted(final IterationEvent e) {
                doTestNormOfResidual(e);
            }

            public void terminationPerformed(final IterationEvent e) {
                doTestNormOfResidual(e);
            }
        };
        solver = new SymmLQ(maxIterations, 1E-10, true);
        solver.getIterationManager().addIterationListener(listener);
        final RealVector b = new ArrayRealVector(n);
        for (int j = 0; j < n; j++) {
            b.set(0.);
            b.setEntry(j, 1.);
            solver.solve(a, b);
        }
    }

[1823, 'src/test/java', 'org.apache.commons.math3.fraction', 'BigFractionTest', 'testReciprocal', 302, 333]

[1823, 'src/test/java', 'org.apache.commons.math3.fraction', 'FractionTest', 'testReciprocal', 274, 304]

    @Test
    public void testReciprocal() {
        BigFraction f = null;

        f = new BigFraction(50, 75);
        f = f.reciprocal();
        Assert.assertEquals(3, f.getNumeratorAsInt());
        Assert.assertEquals(2, f.getDenominatorAsInt());

        f = new BigFraction(4, 3);
        f = f.reciprocal();
        Assert.assertEquals(3, f.getNumeratorAsInt());
        Assert.assertEquals(4, f.getDenominatorAsInt());

        f = new BigFraction(-15, 47);
        f = f.reciprocal();
        Assert.assertEquals(-47, f.getNumeratorAsInt());
        Assert.assertEquals(15, f.getDenominatorAsInt());

        f = new BigFraction(0, 3);
        try {
            f = f.reciprocal();
            Assert.fail("expecting ZeroException");
        } catch (ZeroException ex) {
        }

        // large values
        f = new BigFraction(Integer.MAX_VALUE, 1);
        f = f.reciprocal();
        Assert.assertEquals(1, f.getNumeratorAsInt());
        Assert.assertEquals(Integer.MAX_VALUE, f.getDenominatorAsInt());
    }

    @Test
    public void testReciprocal() {
        Fraction f = null;

        f = new Fraction(50, 75);
        f = f.reciprocal();
        Assert.assertEquals(3, f.getNumerator());
        Assert.assertEquals(2, f.getDenominator());

        f = new Fraction(4, 3);
        f = f.reciprocal();
        Assert.assertEquals(3, f.getNumerator());
        Assert.assertEquals(4, f.getDenominator());

        f = new Fraction(-15, 47);
        f = f.reciprocal();
        Assert.assertEquals(-47, f.getNumerator());
        Assert.assertEquals(15, f.getDenominator());

        f = new Fraction(0, 3);
        try {
            f = f.reciprocal();
            Assert.fail("expecting MathArithmeticException");
        } catch (MathArithmeticException ex) {}

        // large values
        f = new Fraction(Integer.MAX_VALUE, 1);
        f = f.reciprocal();
        Assert.assertEquals(1, f.getNumerator());
        Assert.assertEquals(Integer.MAX_VALUE, f.getDenominator());
    }

[1824, 'src/test/java', 'org.apache.commons.math3.stat.inference', 'TTestTest', 'testTwoSampleTHomoscedastic', 248, 270]

[1824, 'src/test/java', 'org.apache.commons.math3.stat.inference', 'TestUtilsTest', 'testTwoSampleTHomoscedastic', 401, 423]

    @Test
    public void testTwoSampleTHomoscedastic() {
        double[] sample1 ={2, 4, 6, 8, 10, 97};
        double[] sample2 = {4, 6, 8, 10, 16};
        SummaryStatistics sampleStats1 = new SummaryStatistics();
        for (int i = 0; i < sample1.length; i++) {
            sampleStats1.addValue(sample1[i]);
        }
        SummaryStatistics sampleStats2 = new SummaryStatistics();
        for (int i = 0; i < sample2.length; i++) {
            sampleStats2.addValue(sample2[i]);
        }

        // Target comparison values computed using R version 1.8.1 (Linux version)
        Assert.assertEquals("two sample homoscedastic t stat", 0.73096310086,
              testStatistic.homoscedasticT(sample1, sample2), 10E-11);
        Assert.assertEquals("two sample homoscedastic p value", 0.4833963785,
                testStatistic.homoscedasticTTest(sampleStats1, sampleStats2), 1E-10);
        Assert.assertTrue("two sample homoscedastic t-test reject",
                testStatistic.homoscedasticTTest(sample1, sample2, 0.49));
        Assert.assertTrue("two sample homoscedastic t-test accept",
                !testStatistic.homoscedasticTTest(sample1, sample2, 0.48));
    }

    @Test
    public void testTwoSampleTHomoscedastic() {
        double[] sample1 ={2, 4, 6, 8, 10, 97};
        double[] sample2 = {4, 6, 8, 10, 16};
        SummaryStatistics sampleStats1 = new SummaryStatistics();
        for (int i = 0; i < sample1.length; i++) {
            sampleStats1.addValue(sample1[i]);
        }
        SummaryStatistics sampleStats2 = new SummaryStatistics();
        for (int i = 0; i < sample2.length; i++) {
            sampleStats2.addValue(sample2[i]);
        }

        // Target comparison values computed using R version 1.8.1 (Linux version)
        Assert.assertEquals("two sample homoscedastic t stat", 0.73096310086,
                TestUtils.homoscedasticT(sample1, sample2), 10E-11);
        Assert.assertEquals("two sample homoscedastic p value", 0.4833963785,
                TestUtils.homoscedasticTTest(sampleStats1, sampleStats2), 1E-10);
        Assert.assertTrue("two sample homoscedastic t-test reject",
                TestUtils.homoscedasticTTest(sample1, sample2, 0.49));
        Assert.assertTrue("two sample homoscedastic t-test accept",
                !TestUtils.homoscedasticTTest(sample1, sample2, 0.48));
    }

[1832, 'src/test/java', 'org.apache.commons.math3.linear', 'BlockFieldMatrixTest', 'testExamples', 514, 565]

[1832, 'src/test/java', 'org.apache.commons.math3.linear', 'FieldMatrixImplTest', 'testExamples', 377, 425]

    @Test
    public void testExamples() {
        // Create a real matrix with two rows and three columns
        Fraction[][] matrixData = {
                {new Fraction(1),new Fraction(2),new Fraction(3)},
                {new Fraction(2),new Fraction(5),new Fraction(3)}
        };
        FieldMatrix<Fraction> m = new BlockFieldMatrix<Fraction>(matrixData);
        // One more with three rows, two columns
        Fraction[][] matrixData2 = {
                {new Fraction(1),new Fraction(2)},
                {new Fraction(2),new Fraction(5)},
                {new Fraction(1), new Fraction(7)}
        };
        FieldMatrix<Fraction> n = new BlockFieldMatrix<Fraction>(matrixData2);
        // Now multiply m by n
        FieldMatrix<Fraction> p = m.multiply(n);
        Assert.assertEquals(2, p.getRowDimension());
        Assert.assertEquals(2, p.getColumnDimension());
        // Invert p
        FieldMatrix<Fraction> pInverse = new FieldLUDecomposition<Fraction>(p).getSolver().getInverse();
        Assert.assertEquals(2, pInverse.getRowDimension());
        Assert.assertEquals(2, pInverse.getColumnDimension());

        // Solve example
        Fraction[][] coefficientsData = {
                {new Fraction(2), new Fraction(3), new Fraction(-2)},
                {new Fraction(-1), new Fraction(7), new Fraction(6)},
                {new Fraction(4), new Fraction(-3), new Fraction(-5)}
        };
        FieldMatrix<Fraction> coefficients = new BlockFieldMatrix<Fraction>(coefficientsData);
        Fraction[] constants = {
            new Fraction(1), new Fraction(-2), new Fraction(1)
        };
        Fraction[] solution;
        solution = new FieldLUDecomposition<Fraction>(coefficients)
            .getSolver()
            .solve(new ArrayFieldVector<Fraction>(constants, false)).toArray();
        Assert.assertEquals(new Fraction(2).multiply(solution[0]).
                     add(new Fraction(3).multiply(solution[1])).
                     subtract(new Fraction(2).multiply(solution[2])),
                     constants[0]);
        Assert.assertEquals(new Fraction(-1).multiply(solution[0]).
                     add(new Fraction(7).multiply(solution[1])).
                     add(new Fraction(6).multiply(solution[2])),
                     constants[1]);
        Assert.assertEquals(new Fraction(4).multiply(solution[0]).
                     subtract(new Fraction(3).multiply(solution[1])).
                     subtract(new Fraction(5).multiply(solution[2])),
                     constants[2]);

    }

    @Test
    public void testExamples() {
        // Create a real matrix with two rows and three columns
        Fraction[][] matrixData = {
                {new Fraction(1),new Fraction(2),new Fraction(3)},
                {new Fraction(2),new Fraction(5),new Fraction(3)}
        };
        FieldMatrix<Fraction> m = new Array2DRowFieldMatrix<Fraction>(matrixData);
        // One more with three rows, two columns
        Fraction[][] matrixData2 = {
                {new Fraction(1),new Fraction(2)},
                {new Fraction(2),new Fraction(5)},
                {new Fraction(1), new Fraction(7)}
        };
        FieldMatrix<Fraction> n = new Array2DRowFieldMatrix<Fraction>(matrixData2);
        // Now multiply m by n
        FieldMatrix<Fraction> p = m.multiply(n);
        Assert.assertEquals(2, p.getRowDimension());
        Assert.assertEquals(2, p.getColumnDimension());
        // Invert p
        FieldMatrix<Fraction> pInverse = new FieldLUDecomposition<Fraction>(p).getSolver().getInverse();
        Assert.assertEquals(2, pInverse.getRowDimension());
        Assert.assertEquals(2, pInverse.getColumnDimension());

        // Solve example
        Fraction[][] coefficientsData = {
                {new Fraction(2), new Fraction(3), new Fraction(-2)},
                {new Fraction(-1), new Fraction(7), new Fraction(6)},
                {new Fraction(4), new Fraction(-3), new Fraction(-5)}
        };
        FieldMatrix<Fraction> coefficients = new Array2DRowFieldMatrix<Fraction>(coefficientsData);
        Fraction[] constants = {
            new Fraction(1), new Fraction(-2), new Fraction(1)
        };
        Fraction[] solution;
        solution = new FieldLUDecomposition<Fraction>(coefficients)
            .getSolver()
            .solve(new ArrayFieldVector<Fraction>(constants, false)).toArray();
        Assert.assertEquals(new Fraction(2).multiply(solution[0]).
                     add(new Fraction(3).multiply(solution[1])).
                     subtract(new Fraction(2).multiply(solution[2])), constants[0]);
        Assert.assertEquals(new Fraction(-1).multiply(solution[0]).
                     add(new Fraction(7).multiply(solution[1])).
                     add(new Fraction(6).multiply(solution[2])), constants[1]);
        Assert.assertEquals(new Fraction(4).multiply(solution[0]).
                     subtract(new Fraction(3).multiply(solution[1])).
                     subtract(new Fraction(5).multiply(solution[2])), constants[2]);

    }

[1834, 'src/test/java', 'org.apache.commons.math3.stat.inference', 'GTestTest', 'testGTestIndependance1', 94, 112]

[1834, 'src/test/java', 'org.apache.commons.math3.stat.inference', 'GTestTest', 'testGTestIndependance2', 114, 132]

    @Test
    public void testGTestIndependance1() throws Exception {
        final long[] obs1 = new long[]{
            268, 199, 42
        };

        final long[] obs2 = new long[]{
            807, 759, 184
        };

        final double g = testStatistic.gDataSetsComparison(obs1, obs2);

        Assert.assertEquals("G test statistic",
                7.3008170, g, 1E-6);
        final double p_gti = testStatistic.gTestDataSetsComparison(obs1, obs2);

        Assert.assertEquals("g-Test p-value", 0.0259805, p_gti, 1E-6);
        Assert.assertTrue(testStatistic.gTestDataSetsComparison(obs1, obs2, 0.05));
    }

    @Test
    public void testGTestIndependance2() throws Exception {
        final long[] obs1 = new long[]{
            127, 99, 264
        };

        final long[] obs2 = new long[]{
            116, 67, 161
        };

        final double g = testStatistic.gDataSetsComparison(obs1, obs2);

        Assert.assertEquals("G test statistic",
                6.227288, g, 1E-6);
        final double p_gti = testStatistic.gTestDataSetsComparison(obs1, obs2);

        Assert.assertEquals("g-Test p-value", 0.04443, p_gti, 1E-5);
        Assert.assertTrue(testStatistic.gTestDataSetsComparison(obs1, obs2, 0.05));
    }

[1850, 'src/test/java', 'org.apache.commons.math3.linear', 'ConjugateGradientTest', 'testNonSquareOperator', 32, 40]

[1850, 'src/test/java', 'org.apache.commons.math3.linear', 'SymmLQTest', 'testNonSquareOperator', 215, 223]

    @Test(expected = NonSquareOperatorException.class)
    public void testNonSquareOperator() {
        final Array2DRowRealMatrix a = new Array2DRowRealMatrix(2, 3);
        final IterativeLinearSolver solver;
        solver = new ConjugateGradient(10, 0., false);
        final ArrayRealVector b = new ArrayRealVector(a.getRowDimension());
        final ArrayRealVector x = new ArrayRealVector(a.getColumnDimension());
        solver.solve(a, b, x);
    }

    @Test(expected = NonSquareOperatorException.class)
    public void testNonSquareOperator() {
        final Array2DRowRealMatrix a = new Array2DRowRealMatrix(2, 3);
        final IterativeLinearSolver solver;
        solver = new SymmLQ(10, 0., false);
        final ArrayRealVector b = new ArrayRealVector(a.getRowDimension());
        final ArrayRealVector x = new ArrayRealVector(a.getColumnDimension());
        solver.solve(a, b, x);
    }

[1851, 'src/test/java', 'org.apache.commons.math3.linear', 'BlockFieldMatrixTest', 'testPlusMinus', 212, 223]

[1851, 'src/test/java', 'org.apache.commons.math3.linear', 'FieldMatrixImplTest', 'testPlusMinus', 153, 164]

    @Test
    public void testPlusMinus() {
        BlockFieldMatrix<Fraction> m = new BlockFieldMatrix<Fraction>(testData);
        BlockFieldMatrix<Fraction> m2 = new BlockFieldMatrix<Fraction>(testDataInv);
        TestUtils.assertEquals(m.subtract(m2), m2.scalarMultiply(new Fraction(-1)).add(m));
        try {
            m.subtract(new BlockFieldMatrix<Fraction>(testData2));
            Assert.fail("Expecting illegalArgumentException");
        } catch (MathIllegalArgumentException ex) {
            // ignored
        }
    }

    @Test
    public void testPlusMinus() {
        Array2DRowFieldMatrix<Fraction> m = new Array2DRowFieldMatrix<Fraction>(testData);
        Array2DRowFieldMatrix<Fraction> m2 = new Array2DRowFieldMatrix<Fraction>(testDataInv);
        TestUtils.assertEquals(m.subtract(m2),m2.scalarMultiply(new Fraction(-1)).add(m));
        try {
            m.subtract(new Array2DRowFieldMatrix<Fraction>(testData2));
            Assert.fail("Expecting illegalArgumentException");
        } catch (MathIllegalArgumentException ex) {
            // ignored
        }
    }

[1855, 'src/test/java', 'org.apache.commons.math3.distribution', 'ExponentialDistributionTest', 'testCumulativeProbabilityExtremes', 72, 77]

[1855, 'src/test/java', 'org.apache.commons.math3.distribution', 'FDistributionTest', 'testCumulativeProbabilityExtremes', 69, 74]

    @Test
    public void testCumulativeProbabilityExtremes() {
        setCumulativeTestPoints(new double[] {-2, 0});
        setCumulativeTestValues(new double[] {0, 0});
        verifyCumulativeProbabilities();
    }

    @Test
    public void testCumulativeProbabilityExtremes() {
        setCumulativeTestPoints(new double[] {-2, 0});
        setCumulativeTestValues(new double[] {0, 0});
        verifyCumulativeProbabilities();
    }

[1893, 'src/test/java', 'org.apache.commons.math3.linear', 'QRSolverTest', 'testSolve', 112, 143]

[1893, 'src/test/java', 'org.apache.commons.math3.linear', 'RRQRSolverTest', 'testSolve', 114, 146]

    @Test
    public void testSolve() {
        QRDecomposition decomposition =
            new QRDecomposition(MatrixUtils.createRealMatrix(testData3x3NonSingular));
        DecompositionSolver solver = decomposition.getSolver();
        RealMatrix b = MatrixUtils.createRealMatrix(new double[][] {
                { -102, 12250 }, { 544, 24500 }, { 167, -36750 }
        });
        RealMatrix xRef = MatrixUtils.createRealMatrix(new double[][] {
                { 1, 2515 }, { 2, 422 }, { -3, 898 }
        });

        // using RealMatrix
        Assert.assertEquals(0, solver.solve(b).subtract(xRef).getNorm(), 2.0e-16 * xRef.getNorm());

        // using ArrayRealVector
        for (int i = 0; i < b.getColumnDimension(); ++i) {
            final RealVector x = solver.solve(b.getColumnVector(i));
            final double error = x.subtract(xRef.getColumnVector(i)).getNorm();
            Assert.assertEquals(0, error, 3.0e-16 * xRef.getColumnVector(i).getNorm());
        }

        // using RealVector with an alternate implementation
        for (int i = 0; i < b.getColumnDimension(); ++i) {
            ArrayRealVectorTest.RealVectorTestImpl v =
                new ArrayRealVectorTest.RealVectorTestImpl(b.getColumn(i));
            final RealVector x = solver.solve(v);
            final double error = x.subtract(xRef.getColumnVector(i)).getNorm();
            Assert.assertEquals(0, error, 3.0e-16 * xRef.getColumnVector(i).getNorm());
        }

    }

    @Test
    public void testSolve() {
        RealMatrix b = MatrixUtils.createRealMatrix(new double[][] {
                { -102, 12250 }, { 544, 24500 }, { 167, -36750 }
        });
        RealMatrix xRef = MatrixUtils.createRealMatrix(new double[][] {
                { 1, 2515 }, { 2, 422 }, { -3, 898 }
        });


        RRQRDecomposition decomposition = new RRQRDecomposition(MatrixUtils.createRealMatrix(testData3x3NonSingular));
        DecompositionSolver solver = decomposition.getSolver();

        // using RealMatrix
        Assert.assertEquals(0, solver.solve(b).subtract(xRef).getNorm(), 3.0e-16 * xRef.getNorm());

        // using ArrayRealVector
        for (int i = 0; i < b.getColumnDimension(); ++i) {
            final RealVector x = solver.solve(b.getColumnVector(i));
            final double error = x.subtract(xRef.getColumnVector(i)).getNorm();
            Assert.assertEquals(0, error, 3.0e-16 * xRef.getColumnVector(i).getNorm());
        }

        // using RealVector with an alternate implementation
        for (int i = 0; i < b.getColumnDimension(); ++i) {
            ArrayRealVectorTest.RealVectorTestImpl v =
                new ArrayRealVectorTest.RealVectorTestImpl(b.getColumn(i));
            final RealVector x = solver.solve(v);
            final double error = x.subtract(xRef.getColumnVector(i)).getNorm();
            Assert.assertEquals(0, error, 3.0e-16 * xRef.getColumnVector(i).getNorm());
        }

    }

[1898, 'src/test/java', 'org.apache.commons.math3.geometry.euclidean.twod', 'PolygonsSetTest', 'testUnion', 384, 446]

[1898, 'src/test/java', 'org.apache.commons.math3.geometry.euclidean.twod', 'PolygonsSetTest', 'testDifference', 568, 629]

    @Test
    public void testUnion() {
        Vector2D[][] vertices1 = new Vector2D[][] {
            new Vector2D[] {
                new Vector2D( 0.0,  0.0),
                new Vector2D( 2.0,  0.0),
                new Vector2D( 2.0,  2.0),
                new Vector2D( 0.0,  2.0)
            }
        };
        PolygonsSet set1 = buildSet(vertices1);
        Vector2D[][] vertices2 = new Vector2D[][] {
            new Vector2D[] {
                new Vector2D( 1.0,  1.0),
                new Vector2D( 3.0,  1.0),
                new Vector2D( 3.0,  3.0),
                new Vector2D( 1.0,  3.0)
            }
        };
        PolygonsSet set2 = buildSet(vertices2);
        PolygonsSet set  = (PolygonsSet) new RegionFactory<Euclidean2D>().union(set1.copySelf(),
                                                                                set2.copySelf());
        checkVertices(set1.getVertices(), vertices1);
        checkVertices(set2.getVertices(), vertices2);
        checkVertices(set.getVertices(), new Vector2D[][] {
            new Vector2D[] {
                new Vector2D( 0.0,  0.0),
                new Vector2D( 2.0,  0.0),
                new Vector2D( 2.0,  1.0),
                new Vector2D( 3.0,  1.0),
                new Vector2D( 3.0,  3.0),
                new Vector2D( 1.0,  3.0),
                new Vector2D( 1.0,  2.0),
                new Vector2D( 0.0,  2.0)
            }
        });
        checkPoints(Region.Location.INSIDE, set, new Vector2D[] {
            new Vector2D(1.0, 1.0),
            new Vector2D(0.5, 0.5),
            new Vector2D(2.0, 2.0),
            new Vector2D(2.5, 2.5),
            new Vector2D(0.5, 1.5),
            new Vector2D(1.5, 1.5),
            new Vector2D(1.5, 0.5),
            new Vector2D(1.5, 2.5),
            new Vector2D(2.5, 1.5),
            new Vector2D(2.5, 2.5)
        });
        checkPoints(Region.Location.OUTSIDE, set, new Vector2D[] {
            new Vector2D(-0.5, 0.5),
            new Vector2D( 0.5, 2.5),
            new Vector2D( 2.5, 0.5),
            new Vector2D( 3.5, 2.5)
        });
        checkPoints(Region.Location.BOUNDARY, set, new Vector2D[] {
            new Vector2D(0.0, 0.0),
            new Vector2D(0.5, 2.0),
            new Vector2D(2.0, 0.5),
            new Vector2D(2.5, 1.0),
            new Vector2D(3.0, 2.5)
        });

    }

    @Test
    public void testDifference() {
        Vector2D[][] vertices1 = new Vector2D[][] {
            new Vector2D[] {
                new Vector2D( 0.0,  0.0),
                new Vector2D( 2.0,  0.0),
                new Vector2D( 2.0,  2.0),
                new Vector2D( 0.0,  2.0)
            }
        };
        PolygonsSet set1 = buildSet(vertices1);
        Vector2D[][] vertices2 = new Vector2D[][] {
            new Vector2D[] {
                new Vector2D( 1.0,  1.0),
                new Vector2D( 3.0,  1.0),
                new Vector2D( 3.0,  3.0),
                new Vector2D( 1.0,  3.0)
            }
        };
        PolygonsSet set2 = buildSet(vertices2);
        PolygonsSet set  = (PolygonsSet) new RegionFactory<Euclidean2D>().difference(set1.copySelf(),
                                                                                     set2.copySelf());
        checkVertices(set1.getVertices(), vertices1);
        checkVertices(set2.getVertices(), vertices2);
        checkVertices(set.getVertices(), new Vector2D[][] {
            new Vector2D[] {
                new Vector2D( 0.0,  0.0),
                new Vector2D( 2.0,  0.0),
                new Vector2D( 2.0,  1.0),
                new Vector2D( 1.0,  1.0),
                new Vector2D( 1.0,  2.0),
                new Vector2D( 0.0,  2.0)
            }
        });
        checkPoints(Region.Location.INSIDE, set, new Vector2D[] {
            new Vector2D(0.5, 0.5),
            new Vector2D(0.5, 1.5),
            new Vector2D(1.5, 0.5)
        });
        checkPoints(Region.Location.OUTSIDE, set, new Vector2D[] {
            new Vector2D( 2.5, 2.5),
            new Vector2D(-0.5, 0.5),
            new Vector2D( 0.5, 2.5),
            new Vector2D( 2.5, 0.5),
            new Vector2D( 1.5, 1.5),
            new Vector2D( 3.5, 2.5),
            new Vector2D( 1.5, 2.5),
            new Vector2D( 2.5, 1.5),
            new Vector2D( 2.0, 1.5),
            new Vector2D( 2.0, 2.0),
            new Vector2D( 2.5, 1.0),
            new Vector2D( 2.5, 2.5),
            new Vector2D( 3.0, 2.5)
        });
        checkPoints(Region.Location.BOUNDARY, set, new Vector2D[] {
            new Vector2D(1.0, 1.0),
            new Vector2D(1.5, 1.0),
            new Vector2D(0.0, 0.0),
            new Vector2D(0.5, 2.0),
            new Vector2D(2.0, 0.5)
        });
    }

[1901, 'src/test/java', 'org.apache.commons.math3.distribution', 'EnumeratedRealDistributionTest', 'testDensity', 107, 115]

[1901, 'src/test/java', 'org.apache.commons.math3.distribution', 'EnumeratedRealDistributionTest', 'testProbability', 94, 102]

    @Test
    public void testDensity() {
        double[] points = new double[]{-2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0};
        double[] results = new double[]{0, 0.2, 0, 0, 0, 0.5, 0, 0, 0, 0.3, 0};
        for (int p = 0; p < points.length; p++) {
            double density = testDistribution.density(points[p]);
            Assert.assertEquals(results[p], density, 0.0);
        }
    }

    @Test
    public void testProbability() {
        double[] points = new double[]{-2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0};
        double[] results = new double[]{0, 0.2, 0, 0, 0, 0.5, 0, 0, 0, 0.3, 0};
        for (int p = 0; p < points.length; p++) {
            double density = testDistribution.probability(points[p]);
            Assert.assertEquals(results[p], density, 0.0);
        }
    }

[1913, 'src/test/java', 'org.apache.commons.math3.analysis.interpolation', 'BicubicInterpolatorTest', 'testPreconditions', 35, 76]

[1913, 'src/test/java', 'org.apache.commons.math3.analysis.interpolation', 'BicubicSplineInterpolatorTest', 'testPreconditions', 38, 79]

    @Test
    public void testPreconditions() {
        double[] xval = new double[] {3, 4, 5, 6.5};
        double[] yval = new double[] {-4, -3, -1, 2.5};
        double[][] zval = new double[xval.length][yval.length];

        BivariateGridInterpolator interpolator = new BicubicInterpolator();

        @SuppressWarnings("unused")
        BivariateFunction p = interpolator.interpolate(xval, yval, zval);

        double[] wxval = new double[] {3, 2, 5, 6.5};
        try {
            p = interpolator.interpolate(wxval, yval, zval);
            Assert.fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException e) {
            // Expected
        }

        double[] wyval = new double[] {-4, -3, -1, -1};
        try {
            p = interpolator.interpolate(xval, wyval, zval);
            Assert.fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException e) {
            // Expected
        }

        double[][] wzval = new double[xval.length][yval.length + 1];
        try {
            p = interpolator.interpolate(xval, yval, wzval);
            Assert.fail("an exception should have been thrown");
        } catch (DimensionMismatchException e) {
            // Expected
        }
        wzval = new double[xval.length - 1][yval.length];
        try {
            p = interpolator.interpolate(xval, yval, wzval);
            Assert.fail("an exception should have been thrown");
        } catch (DimensionMismatchException e) {
            // Expected
        }
    }

    @Test
    public void testPreconditions() {
        double[] xval = new double[] {3, 4, 5, 6.5};
        double[] yval = new double[] {-4, -3, -1, 2.5};
        double[][] zval = new double[xval.length][yval.length];

        BivariateGridInterpolator interpolator = new BicubicSplineInterpolator();

        @SuppressWarnings("unused")
        BivariateFunction p = interpolator.interpolate(xval, yval, zval);

        double[] wxval = new double[] {3, 2, 5, 6.5};
        try {
            p = interpolator.interpolate(wxval, yval, zval);
            Assert.fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException e) {
            // Expected
        }

        double[] wyval = new double[] {-4, -3, -1, -1};
        try {
            p = interpolator.interpolate(xval, wyval, zval);
            Assert.fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException e) {
            // Expected
        }

        double[][] wzval = new double[xval.length][yval.length + 1];
        try {
            p = interpolator.interpolate(xval, yval, wzval);
            Assert.fail("an exception should have been thrown");
        } catch (DimensionMismatchException e) {
            // Expected
        }
        wzval = new double[xval.length - 1][yval.length];
        try {
            p = interpolator.interpolate(xval, yval, wzval);
            Assert.fail("an exception should have been thrown");
        } catch (DimensionMismatchException e) {
            // Expected
        }
    }

[1932, 'src/test/java', 'org.apache.commons.math3.geometry.euclidean.threed', 'Vector3DTest', 'testNorm1', 134, 138]

[1932, 'src/test/java', 'org.apache.commons.math3.geometry.euclidean.threed', 'Vector3DTest', 'testNormInf', 152, 156]

    @Test
    public void testNorm1() {
        Assert.assertEquals(0.0, Vector3D.ZERO.getNorm1(), 0);
        Assert.assertEquals(6.0, new Vector3D(1, -2, 3).getNorm1(), 0);
    }

    @Test
    public void testNormInf() {
        Assert.assertEquals(0.0, Vector3D.ZERO.getNormInf(), 0);
        Assert.assertEquals(3.0, new Vector3D(1, -2, 3).getNormInf(), 0);
    }

[1933, 'src/test/java', 'org.apache.commons.math3.geometry.euclidean.oned', 'Vector1DFormatAbstractTest', 'testParseNegativeX', 170, 178]

[1933, 'src/test/java', 'org.apache.commons.math3.geometry.euclidean.oned', 'Vector1DFormatAbstractTest', 'testParseNegativeAll', 200, 208]

    @Test
    public void testParseNegativeX() throws MathParseException {
        String source =
            "{-1" + getDecimalCharacter() +
            "2323}";
        Vector1D expected = new Vector1D(-1.2323);
        Vector1D actual = vector1DFormat.parse(source);
        Assert.assertEquals(expected, actual);
    }

    @Test
    public void testParseNegativeAll() throws MathParseException {
        String source =
            "{-1" + getDecimalCharacter() +
            "2323}";
        Vector1D expected = new Vector1D(-1.2323);
        Vector1D actual = vector1DFormat.parse(source);
        Assert.assertEquals(expected, actual);
    }

[1935, 'src/test/java', 'org.apache.commons.math3.linear', 'ConjugateGradientTest', 'testUnpreconditionedSolution', 78, 99]

[1935, 'src/test/java', 'org.apache.commons.math3.linear', 'SymmLQTest', 'testUnpreconditionedSolution', 244, 265]

    @Test
    public void testUnpreconditionedSolution() {
        final int n = 5;
        final int maxIterations = 100;
        final RealLinearOperator a = new HilbertMatrix(n);
        final InverseHilbertMatrix ainv = new InverseHilbertMatrix(n);
        final IterativeLinearSolver solver;
        solver = new ConjugateGradient(maxIterations, 1E-10, true);
        final RealVector b = new ArrayRealVector(n);
        for (int j = 0; j < n; j++) {
            b.set(0.);
            b.setEntry(j, 1.);
            final RealVector x = solver.solve(a, b);
            for (int i = 0; i < n; i++) {
                final double actual = x.getEntry(i);
                final double expected = ainv.getEntry(i, j);
                final double delta = 1E-10 * FastMath.abs(expected);
                final String msg = String.format("entry[%d][%d]", i, j);
                Assert.assertEquals(msg, expected, actual, delta);
            }
        }
    }

    @Test
    public void testUnpreconditionedSolution() {
        final int n = 5;
        final int maxIterations = 100;
        final RealLinearOperator a = new HilbertMatrix(n);
        final InverseHilbertMatrix ainv = new InverseHilbertMatrix(n);
        final IterativeLinearSolver solver;
        solver = new SymmLQ(maxIterations, 1E-10, true);
        final RealVector b = new ArrayRealVector(n);
        for (int j = 0; j < n; j++) {
            b.set(0.);
            b.setEntry(j, 1.);
            final RealVector x = solver.solve(a, b);
            for (int i = 0; i < n; i++) {
                final double actual = x.getEntry(i);
                final double expected = ainv.getEntry(i, j);
                final double delta = 1E-6 * FastMath.abs(expected);
                final String msg = String.format("entry[%d][%d]", i, j);
                Assert.assertEquals(msg, expected, actual, delta);
            }
        }
    }

[1955, 'src/test/java', 'org.apache.commons.math3.transform', 'FastFourierTransformerTest', 'testTransformFunctionSizeNotAPowerOfTwo', 93, 114]

[1955, 'src/test/java', 'org.apache.commons.math3.transform', 'FastFourierTransformerTest', 'testTransformFunctionInvalidBounds', 140, 161]

    @Test
    public void testTransformFunctionSizeNotAPowerOfTwo() {
        final int n = 127;
        final UnivariateFunction f = new Sin();
        final DftNormalization[] norm;
        norm = DftNormalization.values();
        final TransformType[] type;
        type = TransformType.values();
        for (int i = 0; i < norm.length; i++) {
            for (int j = 0; j < type.length; j++) {
                final FastFourierTransformer fft;
                fft = new FastFourierTransformer(norm[i]);
                try {
                    fft.transform(f, 0.0, Math.PI, n, type[j]);
                    Assert.fail(norm[i] + ", " + type[j] +
                        ": MathIllegalArgumentException was expected");
                } catch (MathIllegalArgumentException e) {
                    // Expected behaviour
                }
            }
        }
    }

    @Test
    public void testTransformFunctionInvalidBounds() {
        final int n = 128;
        final UnivariateFunction f = new Sin();
        final DftNormalization[] norm;
        norm = DftNormalization.values();
        final TransformType[] type;
        type = TransformType.values();
        for (int i = 0; i < norm.length; i++) {
            for (int j = 0; j < type.length; j++) {
                final FastFourierTransformer fft;
                fft = new FastFourierTransformer(norm[i]);
                try {
                    fft.transform(f, Math.PI, 0.0, n, type[j]);
                    Assert.fail(norm[i] + ", " + type[j] +
                        ": NumberIsTooLargeException was expected");
                } catch (NumberIsTooLargeException e) {
                    // Expected behaviour
                }
            }
        }
    }

[1956, 'src/test/java', 'org.apache.commons.math3.stat.descriptive.rank', 'MaxTest', 'testNaNs', 64, 72]

[1956, 'src/test/java', 'org.apache.commons.math3.stat.descriptive.rank', 'MinTest', 'testNaNs', 64, 72]

    @Test
    public void testNaNs() {
        Max max = new Max();
        double nan = Double.NaN;
        Assert.assertEquals(3d, max.evaluate(new double[]{nan, 2d, 3d}), 0);
        Assert.assertEquals(3d, max.evaluate(new double[]{1d, nan, 3d}), 0);
        Assert.assertEquals(2d, max.evaluate(new double[]{1d, 2d, nan}), 0);
        Assert.assertTrue(Double.isNaN(max.evaluate(new double[]{nan, nan, nan})));
    }

    @Test
    public void testNaNs() {
        Min min = new Min();
        double nan = Double.NaN;
        Assert.assertEquals(2d, min.evaluate(new double[]{nan, 2d, 3d}), 0);
        Assert.assertEquals(1d, min.evaluate(new double[]{1d, nan, 3d}), 0);
        Assert.assertEquals(1d, min.evaluate(new double[]{1d, 2d, nan}), 0);
        Assert.assertTrue(Double.isNaN(min.evaluate(new double[]{nan, nan, nan})));
    }

[1963, 'src/test/java', 'org.apache.commons.math3.complex', 'ComplexFormatAbstractTest', 'testGetImaginaryFormat', 282, 287]

[1963, 'src/test/java', 'org.apache.commons.math3.complex', 'ComplexFormatAbstractTest', 'testGetRealFormat', 289, 294]

    @Test
    public void testGetImaginaryFormat() {
        NumberFormat nf = NumberFormat.getInstance();
        ComplexFormat cf = new ComplexFormat(nf);
        Assert.assertSame(nf, cf.getImaginaryFormat());
    }

    @Test
    public void testGetRealFormat() {
        NumberFormat nf = NumberFormat.getInstance();
        ComplexFormat cf = new ComplexFormat(nf);
        Assert.assertSame(nf, cf.getRealFormat());
    }

[1969, 'src/test/java', 'org.apache.commons.math3.fraction', 'BigFractionFormatTest', 'testParseNegative', 175, 201]

[1969, 'src/test/java', 'org.apache.commons.math3.fraction', 'FractionFormatTest', 'testParseNegative', 213, 239]

    @Test
    public void testParseNegative() {

        {
            String source = "-1 / 2";
            BigFraction c = properFormat.parse(source);
            Assert.assertNotNull(c);
            Assert.assertEquals(-1, c.getNumeratorAsInt());
            Assert.assertEquals(2, c.getDenominatorAsInt());

            c = improperFormat.parse(source);
            Assert.assertNotNull(c);
            Assert.assertEquals(-1, c.getNumeratorAsInt());
            Assert.assertEquals(2, c.getDenominatorAsInt());

            source = "1 / -2";
            c = properFormat.parse(source);
            Assert.assertNotNull(c);
            Assert.assertEquals(-1, c.getNumeratorAsInt());
            Assert.assertEquals(2, c.getDenominatorAsInt());

            c = improperFormat.parse(source);
            Assert.assertNotNull(c);
            Assert.assertEquals(-1, c.getNumeratorAsInt());
            Assert.assertEquals(2, c.getDenominatorAsInt());
        }
    }

    @Test
    public void testParseNegative() {

        {
            String source = "-1 / 2";
            Fraction c = properFormat.parse(source);
            Assert.assertNotNull(c);
            Assert.assertEquals(-1, c.getNumerator());
            Assert.assertEquals(2, c.getDenominator());

            c = improperFormat.parse(source);
            Assert.assertNotNull(c);
            Assert.assertEquals(-1, c.getNumerator());
            Assert.assertEquals(2, c.getDenominator());

            source = "1 / -2";
            c = properFormat.parse(source);
            Assert.assertNotNull(c);
            Assert.assertEquals(-1, c.getNumerator());
            Assert.assertEquals(2, c.getDenominator());

            c = improperFormat.parse(source);
            Assert.assertNotNull(c);
            Assert.assertEquals(-1, c.getNumerator());
            Assert.assertEquals(2, c.getDenominator());
        }
    }

[1987, 'src/test/java', 'org.apache.commons.math3.linear', 'BlockRealMatrixTest', 'testOperateLarge', 340, 352]

[1987, 'src/test/java', 'org.apache.commons.math3.linear', 'BlockRealMatrixTest', 'testOperatePremultiplyLarge', 354, 366]

    @Test
    public void testOperateLarge() {
        int p = (7 * BlockRealMatrix.BLOCK_SIZE) / 2;
        int q = (5 * BlockRealMatrix.BLOCK_SIZE) / 2;
        int r =  3 * BlockRealMatrix.BLOCK_SIZE;
        Random random = new Random(111007463902334l);
        RealMatrix m1 = createRandomMatrix(random, p, q);
        RealMatrix m2 = createRandomMatrix(random, q, r);
        RealMatrix m1m2 = m1.multiply(m2);
        for (int i = 0; i < r; ++i) {
            checkArrays(m1m2.getColumn(i), m1.operate(m2.getColumn(i)));
        }
    }

    @Test
    public void testOperatePremultiplyLarge() {
        int p = (7 * BlockRealMatrix.BLOCK_SIZE) / 2;
        int q = (5 * BlockRealMatrix.BLOCK_SIZE) / 2;
        int r =  3 * BlockRealMatrix.BLOCK_SIZE;
        Random random = new Random(111007463902334l);
        RealMatrix m1 = createRandomMatrix(random, p, q);
        RealMatrix m2 = createRandomMatrix(random, q, r);
        RealMatrix m1m2 = m1.multiply(m2);
        for (int i = 0; i < p; ++i) {
            checkArrays(m1m2.getRow(i), m2.preMultiply(m1.getRow(i)));
        }
    }

[2009, 'src/test/java', 'org.apache.commons.math3.fraction', 'BigFractionTest', 'testCompareTo', 186, 206]

[2009, 'src/test/java', 'org.apache.commons.math3.fraction', 'FractionTest', 'testCompareTo', 167, 186]

    @Test
    public void testCompareTo() {
        BigFraction first = new BigFraction(1, 2);
        BigFraction second = new BigFraction(1, 3);
        BigFraction third = new BigFraction(1, 2);

        Assert.assertEquals(0, first.compareTo(first));
        Assert.assertEquals(0, first.compareTo(third));
        Assert.assertEquals(1, first.compareTo(second));
        Assert.assertEquals(-1, second.compareTo(first));

        // these two values are different approximations of PI
        // the first  one is approximately PI - 3.07e-18
        // the second one is approximately PI + 1.936e-17
        BigFraction pi1 = new BigFraction(1068966896, 340262731);
        BigFraction pi2 = new BigFraction( 411557987, 131002976);
        Assert.assertEquals(-1, pi1.compareTo(pi2));
        Assert.assertEquals( 1, pi2.compareTo(pi1));
        Assert.assertEquals(0.0, pi1.doubleValue() - pi2.doubleValue(), 1.0e-20);

    }

    @Test
    public void testCompareTo() {
        Fraction first = new Fraction(1, 2);
        Fraction second = new Fraction(1, 3);
        Fraction third = new Fraction(1, 2);

        Assert.assertEquals(0, first.compareTo(first));
        Assert.assertEquals(0, first.compareTo(third));
        Assert.assertEquals(1, first.compareTo(second));
        Assert.assertEquals(-1, second.compareTo(first));

        // these two values are different approximations of PI
        // the first  one is approximately PI - 3.07e-18
        // the second one is approximately PI + 1.936e-17
        Fraction pi1 = new Fraction(1068966896, 340262731);
        Fraction pi2 = new Fraction( 411557987, 131002976);
        Assert.assertEquals(-1, pi1.compareTo(pi2));
        Assert.assertEquals( 1, pi2.compareTo(pi1));
        Assert.assertEquals(0.0, pi1.doubleValue() - pi2.doubleValue(), 1.0e-20);
    }

[2023, 'src/test/java', 'org.apache.commons.math3.fraction', 'BigFractionTest', 'testEpsilonLimitConstructor', 174, 184]

[2023, 'src/test/java', 'org.apache.commons.math3.fraction', 'FractionTest', 'testEpsilonLimitConstructor', 155, 165]

    @Test
    public void testEpsilonLimitConstructor() throws ConvergenceException {
        assertFraction(2, 5, new BigFraction(0.4, 1.0e-5, 100));

        assertFraction(3, 5, new BigFraction(0.6152, 0.02, 100));
        assertFraction(8, 13, new BigFraction(0.6152, 1.0e-3, 100));
        assertFraction(251, 408, new BigFraction(0.6152, 1.0e-4, 100));
        assertFraction(251, 408, new BigFraction(0.6152, 1.0e-5, 100));
        assertFraction(510, 829, new BigFraction(0.6152, 1.0e-6, 100));
        assertFraction(769, 1250, new BigFraction(0.6152, 1.0e-7, 100));
    }

    @Test
    public void testEpsilonLimitConstructor() throws ConvergenceException  {
        assertFraction(2, 5, new Fraction(0.4, 1.0e-5, 100));

        assertFraction(3, 5,      new Fraction(0.6152, 0.02, 100));
        assertFraction(8, 13,     new Fraction(0.6152, 1.0e-3, 100));
        assertFraction(251, 408,  new Fraction(0.6152, 1.0e-4, 100));
        assertFraction(251, 408,  new Fraction(0.6152, 1.0e-5, 100));
        assertFraction(510, 829,  new Fraction(0.6152, 1.0e-6, 100));
        assertFraction(769, 1250, new Fraction(0.6152, 1.0e-7, 100));
    }

[2024, 'src/test/java', 'org.apache.commons.math3.random', 'HaltonSequenceGeneratorTest', 'testConstructor', 81, 96]

[2024, 'src/test/java', 'org.apache.commons.math3.random', 'SobolSequenceGeneratorTest', 'testConstructor', 58, 73]

    @Test
    public void testConstructor() {
        try {
            new HaltonSequenceGenerator(0);
            Assert.fail("an exception should have been thrown");
        } catch (OutOfRangeException e) {
            // expected
        }

        try {
            new HaltonSequenceGenerator(41);
            Assert.fail("an exception should have been thrown");
        } catch (OutOfRangeException e) {
            // expected
        }
    }

    @Test
    public void testConstructor() {
        try {
            new SobolSequenceGenerator(0);
            Assert.fail("an exception should have been thrown");
        } catch (OutOfRangeException e) {
            // expected
        }

        try {
            new SobolSequenceGenerator(1001);
            Assert.fail("an exception should have been thrown");
        } catch (OutOfRangeException e) {
            // expected
        }
    }

[2070, 'src/test/java', 'org.apache.commons.math3.linear', 'Array2DRowRealMatrixTest', 'testWalk', 991, 1075]

[2070, 'src/test/java', 'org.apache.commons.math3.linear', 'BlockRealMatrixTest', 'testWalk', 1142, 1227]

    @Test
    public void testWalk() {
        int rows    = 150;
        int columns = 75;

        RealMatrix m = new Array2DRowRealMatrix(rows, columns);
        m.walkInRowOrder(new SetVisitor());
        GetVisitor getVisitor = new GetVisitor();
        m.walkInOptimizedOrder(getVisitor);
        Assert.assertEquals(rows * columns, getVisitor.getCount());

        m = new Array2DRowRealMatrix(rows, columns);
        m.walkInRowOrder(new SetVisitor(), 1, rows - 2, 1, columns - 2);
        getVisitor = new GetVisitor();
        m.walkInOptimizedOrder(getVisitor, 1, rows - 2, 1, columns - 2);
        Assert.assertEquals((rows - 2) * (columns - 2), getVisitor.getCount());
        for (int i = 0; i < rows; ++i) {
            Assert.assertEquals(0.0, m.getEntry(i, 0), 0);
            Assert.assertEquals(0.0, m.getEntry(i, columns - 1), 0);
        }
        for (int j = 0; j < columns; ++j) {
            Assert.assertEquals(0.0, m.getEntry(0, j), 0);
            Assert.assertEquals(0.0, m.getEntry(rows - 1, j), 0);
        }

        m = new Array2DRowRealMatrix(rows, columns);
        m.walkInColumnOrder(new SetVisitor());
        getVisitor = new GetVisitor();
        m.walkInOptimizedOrder(getVisitor);
        Assert.assertEquals(rows * columns, getVisitor.getCount());

        m = new Array2DRowRealMatrix(rows, columns);
        m.walkInColumnOrder(new SetVisitor(), 1, rows - 2, 1, columns - 2);
        getVisitor = new GetVisitor();
        m.walkInOptimizedOrder(getVisitor, 1, rows - 2, 1, columns - 2);
        Assert.assertEquals((rows - 2) * (columns - 2), getVisitor.getCount());
        for (int i = 0; i < rows; ++i) {
            Assert.assertEquals(0.0, m.getEntry(i, 0), 0);
            Assert.assertEquals(0.0, m.getEntry(i, columns - 1), 0);
        }
        for (int j = 0; j < columns; ++j) {
            Assert.assertEquals(0.0, m.getEntry(0, j), 0);
            Assert.assertEquals(0.0, m.getEntry(rows - 1, j), 0);
        }

        m = new Array2DRowRealMatrix(rows, columns);
        m.walkInOptimizedOrder(new SetVisitor());
        getVisitor = new GetVisitor();
        m.walkInRowOrder(getVisitor);
        Assert.assertEquals(rows * columns, getVisitor.getCount());

        m = new Array2DRowRealMatrix(rows, columns);
        m.walkInOptimizedOrder(new SetVisitor(), 1, rows - 2, 1, columns - 2);
        getVisitor = new GetVisitor();
        m.walkInRowOrder(getVisitor, 1, rows - 2, 1, columns - 2);
        Assert.assertEquals((rows - 2) * (columns - 2), getVisitor.getCount());
        for (int i = 0; i < rows; ++i) {
            Assert.assertEquals(0.0, m.getEntry(i, 0), 0);
            Assert.assertEquals(0.0, m.getEntry(i, columns - 1), 0);
        }
        for (int j = 0; j < columns; ++j) {
            Assert.assertEquals(0.0, m.getEntry(0, j), 0);
            Assert.assertEquals(0.0, m.getEntry(rows - 1, j), 0);
        }

        m = new Array2DRowRealMatrix(rows, columns);
        m.walkInOptimizedOrder(new SetVisitor());
        getVisitor = new GetVisitor();
        m.walkInColumnOrder(getVisitor);
        Assert.assertEquals(rows * columns, getVisitor.getCount());

        m = new Array2DRowRealMatrix(rows, columns);
        m.walkInOptimizedOrder(new SetVisitor(), 1, rows - 2, 1, columns - 2);
        getVisitor = new GetVisitor();
        m.walkInColumnOrder(getVisitor, 1, rows - 2, 1, columns - 2);
        Assert.assertEquals((rows - 2) * (columns - 2), getVisitor.getCount());
        for (int i = 0; i < rows; ++i) {
            Assert.assertEquals(0.0, m.getEntry(i, 0), 0);
            Assert.assertEquals(0.0, m.getEntry(i, columns - 1), 0);
        }
        for (int j = 0; j < columns; ++j) {
            Assert.assertEquals(0.0, m.getEntry(0, j), 0);
            Assert.assertEquals(0.0, m.getEntry(rows - 1, j), 0);
        }
    }

    @Test
    public void testWalk() {
        int rows    = 150;
        int columns = 75;

        RealMatrix m = new BlockRealMatrix(rows, columns);
        m.walkInRowOrder(new SetVisitor());
        GetVisitor getVisitor = new GetVisitor();
        m.walkInOptimizedOrder(getVisitor);
        Assert.assertEquals(rows * columns, getVisitor.getCount());

        m = new BlockRealMatrix(rows, columns);
        m.walkInRowOrder(new SetVisitor(), 1, rows - 2, 1, columns - 2);
        getVisitor = new GetVisitor();
        m.walkInOptimizedOrder(getVisitor, 1, rows - 2, 1, columns - 2);
        Assert.assertEquals((rows - 2) * (columns - 2), getVisitor.getCount());
        for (int i = 0; i < rows; ++i) {
            Assert.assertEquals(0.0, m.getEntry(i, 0), 0);
            Assert.assertEquals(0.0, m.getEntry(i, columns - 1), 0);
        }
        for (int j = 0; j < columns; ++j) {
            Assert.assertEquals(0.0, m.getEntry(0, j), 0);
            Assert.assertEquals(0.0, m.getEntry(rows - 1, j), 0);
        }

        m = new BlockRealMatrix(rows, columns);
        m.walkInColumnOrder(new SetVisitor());
        getVisitor = new GetVisitor();
        m.walkInOptimizedOrder(getVisitor);
        Assert.assertEquals(rows * columns, getVisitor.getCount());

        m = new BlockRealMatrix(rows, columns);
        m.walkInColumnOrder(new SetVisitor(), 1, rows - 2, 1, columns - 2);
        getVisitor = new GetVisitor();
        m.walkInOptimizedOrder(getVisitor, 1, rows - 2, 1, columns - 2);
        Assert.assertEquals((rows - 2) * (columns - 2), getVisitor.getCount());
        for (int i = 0; i < rows; ++i) {
            Assert.assertEquals(0.0, m.getEntry(i, 0), 0);
            Assert.assertEquals(0.0, m.getEntry(i, columns - 1), 0);
        }
        for (int j = 0; j < columns; ++j) {
            Assert.assertEquals(0.0, m.getEntry(0, j), 0);
            Assert.assertEquals(0.0, m.getEntry(rows - 1, j), 0);
        }

        m = new BlockRealMatrix(rows, columns);
        m.walkInOptimizedOrder(new SetVisitor());
        getVisitor = new GetVisitor();
        m.walkInRowOrder(getVisitor);
        Assert.assertEquals(rows * columns, getVisitor.getCount());

        m = new BlockRealMatrix(rows, columns);
        m.walkInOptimizedOrder(new SetVisitor(), 1, rows - 2, 1, columns - 2);
        getVisitor = new GetVisitor();
        m.walkInRowOrder(getVisitor, 1, rows - 2, 1, columns - 2);
        Assert.assertEquals((rows - 2) * (columns - 2), getVisitor.getCount());
        for (int i = 0; i < rows; ++i) {
            Assert.assertEquals(0.0, m.getEntry(i, 0), 0);
            Assert.assertEquals(0.0, m.getEntry(i, columns - 1), 0);
        }
        for (int j = 0; j < columns; ++j) {
            Assert.assertEquals(0.0, m.getEntry(0, j), 0);
            Assert.assertEquals(0.0, m.getEntry(rows - 1, j), 0);
        }

        m = new BlockRealMatrix(rows, columns);
        m.walkInOptimizedOrder(new SetVisitor());
        getVisitor = new GetVisitor();
        m.walkInColumnOrder(getVisitor);
        Assert.assertEquals(rows * columns, getVisitor.getCount());

        m = new BlockRealMatrix(rows, columns);
        m.walkInOptimizedOrder(new SetVisitor(), 1, rows - 2, 1, columns - 2);
        getVisitor = new GetVisitor();
        m.walkInColumnOrder(getVisitor, 1, rows - 2, 1, columns - 2);
        Assert.assertEquals((rows - 2) * (columns - 2), getVisitor.getCount());
        for (int i = 0; i < rows; ++i) {
            Assert.assertEquals(0.0, m.getEntry(i, 0), 0);
            Assert.assertEquals(0.0, m.getEntry(i, columns - 1), 0);
        }
        for (int j = 0; j < columns; ++j) {
            Assert.assertEquals(0.0, m.getEntry(0, j), 0);
            Assert.assertEquals(0.0, m.getEntry(rows - 1, j), 0);
        }

    }

[2088, 'src/test/java', 'org.apache.commons.math3.linear', 'BlockFieldMatrixTest', 'testSetRow', 1020, 1038]

[2088, 'src/test/java', 'org.apache.commons.math3.linear', 'FieldMatrixImplTest', 'testSetRow', 780, 798]

    @Test
    public void testSetRow() {
        FieldMatrix<Fraction> m = new BlockFieldMatrix<Fraction>(subTestData);
        Assert.assertTrue(subRow3[0][0] != m.getRow(0)[0]);
        m.setRow(0, subRow3[0]);
        checkArrays(subRow3[0], m.getRow(0));
        try {
            m.setRow(-1, subRow3[0]);
            Assert.fail("Expecting OutOfRangeException");
        } catch (OutOfRangeException ex) {
            // expected
        }
        try {
            m.setRow(0, new Fraction[5]);
            Assert.fail("Expecting MatrixDimensionMismatchException");
        } catch (MatrixDimensionMismatchException ex) {
            // expected
        }
    }

    @Test
    public void testSetRow() {
        FieldMatrix<Fraction> m = new Array2DRowFieldMatrix<Fraction>(subTestData);
        Assert.assertTrue(subRow3[0][0] != m.getRow(0)[0]);
        m.setRow(0, subRow3[0]);
        checkArrays(subRow3[0], m.getRow(0));
        try {
            m.setRow(-1, subRow3[0]);
            Assert.fail("Expecting OutOfRangeException");
        } catch (OutOfRangeException ex) {
            // expected
        }
        try {
            m.setRow(0, new Fraction[5]);
            Assert.fail("Expecting MatrixDimensionMismatchException");
        } catch (MatrixDimensionMismatchException ex) {
            // expected
        }
    }

[2103, 'src/test/java', 'org.apache.commons.math3.analysis.differentiation', 'DerivativeStructureTest', 'testPrimitiveAdd', 101, 108]

[2103, 'src/test/java', 'org.apache.commons.math3.analysis.differentiation', 'DerivativeStructureTest', 'testPrimitiveMultiply', 141, 148]

    @Test
    public void testPrimitiveAdd() {
        for (int maxOrder = 1; maxOrder < 5; ++maxOrder) {
            checkF0F1(new DerivativeStructure(3, maxOrder, 0, 1.0).add(5), 6.0, 1.0, 0.0, 0.0);
            checkF0F1(new DerivativeStructure(3, maxOrder, 1, 2.0).add(5), 7.0, 0.0, 1.0, 0.0);
            checkF0F1(new DerivativeStructure(3, maxOrder, 2, 3.0).add(5), 8.0, 0.0, 0.0, 1.0);
        }
    }

    @Test
    public void testPrimitiveMultiply() {
        for (int maxOrder = 1; maxOrder < 5; ++maxOrder) {
            checkF0F1(new DerivativeStructure(3, maxOrder, 0, 1.0).multiply(5),  5.0, 5.0, 0.0, 0.0);
            checkF0F1(new DerivativeStructure(3, maxOrder, 1, 2.0).multiply(5), 10.0, 0.0, 5.0, 0.0);
            checkF0F1(new DerivativeStructure(3, maxOrder, 2, 3.0).multiply(5), 15.0, 0.0, 0.0, 5.0);
        }
    }

[2112, 'src/test/java', 'org.apache.commons.math3.ode.nonstiff', 'GillIntegratorTest', 'testDecreasingSteps', 55, 101]

[2112, 'src/test/java', 'org.apache.commons.math3.ode.nonstiff', 'ThreeEighthesIntegratorTest', 'testDecreasingSteps', 55, 101]

  @Test
  public void testDecreasingSteps()
      throws DimensionMismatchException, NumberIsTooSmallException,
             MaxCountExceededException, NoBracketingException {

      for (TestProblemAbstract pb : new TestProblemAbstract[] {
          new TestProblem1(), new TestProblem2(), new TestProblem3(),
          new TestProblem4(), new TestProblem5(), new TestProblem6()
      }) {

      double previousValueError = Double.NaN;
      double previousTimeError = Double.NaN;
      for (int i = 5; i < 10; ++i) {

        double step = (pb.getFinalTime() - pb.getInitialTime()) * FastMath.pow(2.0, -i);

        FirstOrderIntegrator integ = new GillIntegrator(step);
        TestProblemHandler handler = new TestProblemHandler(pb, integ);
        integ.addStepHandler(handler);
        EventHandler[] functions = pb.getEventsHandlers();
        for (int l = 0; l < functions.length; ++l) {
          integ.addEventHandler(functions[l],
                                     Double.POSITIVE_INFINITY, 1.0e-6 * step, 1000);
        }
        double stopTime = integ.integrate(pb, pb.getInitialTime(), pb.getInitialState(),
                                          pb.getFinalTime(), new double[pb.getDimension()]);
        if (functions.length == 0) {
            Assert.assertEquals(pb.getFinalTime(), stopTime, 1.0e-10);
        }

        double valueError = handler.getMaximalValueError();
        if (i > 5) {
          Assert.assertTrue(valueError < 1.01 * FastMath.abs(previousValueError));
        }
        previousValueError = valueError;

        double timeError = handler.getMaximalTimeError();
        if (i > 5) {
          Assert.assertTrue(timeError <= FastMath.abs(previousTimeError));
        }
        previousTimeError = timeError;

      }

    }

  }

  @Test
  public void testDecreasingSteps()
      throws DimensionMismatchException, NumberIsTooSmallException,
             MaxCountExceededException, NoBracketingException {

      for (TestProblemAbstract pb : new TestProblemAbstract[] {
          new TestProblem1(), new TestProblem2(), new TestProblem3(),
          new TestProblem4(), new TestProblem5(), new TestProblem6()
      }) {

      double previousValueError = Double.NaN;
      double previousTimeError = Double.NaN;
      for (int i = 4; i < 10; ++i) {

        double step = (pb.getFinalTime() - pb.getInitialTime()) * FastMath.pow(2.0, -i);

        FirstOrderIntegrator integ = new ThreeEighthesIntegrator(step);
        TestProblemHandler handler = new TestProblemHandler(pb, integ);
        integ.addStepHandler(handler);
        EventHandler[] functions = pb.getEventsHandlers();
        for (int l = 0; l < functions.length; ++l) {
          integ.addEventHandler(functions[l],
                                     Double.POSITIVE_INFINITY, 1.0e-6 * step, 1000);
        }
        double stopTime = integ.integrate(pb, pb.getInitialTime(), pb.getInitialState(),
                                          pb.getFinalTime(), new double[pb.getDimension()]);
        if (functions.length == 0) {
            Assert.assertEquals(pb.getFinalTime(), stopTime, 1.0e-10);
        }

        double error = handler.getMaximalValueError();
        if (i > 4) {
          Assert.assertTrue(error < 1.01 * FastMath.abs(previousValueError));
        }
        previousValueError = error;

        double timeError = handler.getMaximalTimeError();
        if (i > 4) {
          Assert.assertTrue(timeError <= FastMath.abs(previousTimeError));
        }
        previousTimeError = timeError;

      }

    }

  }

[2113, 'src/test/java', 'org.apache.commons.math3.linear', 'ConjugateGradientTest', 'testUnpreconditionedSolutionWithInitialGuess', 127, 152]

[2113, 'src/test/java', 'org.apache.commons.math3.linear', 'SymmLQTest', 'testUnpreconditionedSolutionWithInitialGuess', 293, 318]

    @Test
    public void testUnpreconditionedSolutionWithInitialGuess() {
        final int n = 5;
        final int maxIterations = 100;
        final RealLinearOperator a = new HilbertMatrix(n);
        final InverseHilbertMatrix ainv = new InverseHilbertMatrix(n);
        final IterativeLinearSolver solver;
        solver = new ConjugateGradient(maxIterations, 1E-10, true);
        final RealVector b = new ArrayRealVector(n);
        for (int j = 0; j < n; j++) {
            b.set(0.);
            b.setEntry(j, 1.);
            final RealVector x0 = new ArrayRealVector(n);
            x0.set(1.);
            final RealVector x = solver.solve(a, b, x0);
            Assert.assertNotSame("x should not be a reference to x0", x0, x);
            for (int i = 0; i < n; i++) {
                final double actual = x.getEntry(i);
                final double expected = ainv.getEntry(i, j);
                final double delta = 1E-10 * FastMath.abs(expected);
                final String msg = String.format("entry[%d][%d]", i, j);
                Assert.assertEquals(msg, expected, actual, delta);
                Assert.assertEquals(msg, x0.getEntry(i), 1., Math.ulp(1.));
            }
        }
    }

    @Test
    public void testUnpreconditionedSolutionWithInitialGuess() {
        final int n = 5;
        final int maxIterations = 100;
        final RealLinearOperator a = new HilbertMatrix(n);
        final InverseHilbertMatrix ainv = new InverseHilbertMatrix(n);
        final IterativeLinearSolver solver;
        solver = new SymmLQ(maxIterations, 1E-10, true);
        final RealVector b = new ArrayRealVector(n);
        for (int j = 0; j < n; j++) {
            b.set(0.);
            b.setEntry(j, 1.);
            final RealVector x0 = new ArrayRealVector(n);
            x0.set(1.);
            final RealVector x = solver.solve(a, b, x0);
            Assert.assertNotSame("x should not be a reference to x0", x0, x);
            for (int i = 0; i < n; i++) {
                final double actual = x.getEntry(i);
                final double expected = ainv.getEntry(i, j);
                final double delta = 1E-6 * FastMath.abs(expected);
                final String msg = String.format("entry[%d][%d]", i, j);
                Assert.assertEquals(msg, expected, actual, delta);
                Assert.assertEquals(msg, x0.getEntry(i), 1., Math.ulp(1.));
            }
        }
    }

[2124, 'src/test/java', 'org.apache.commons.math3.complex', 'ComplexTest', 'testDivideZero', 226, 232]

[2124, 'src/test/java', 'org.apache.commons.math3.complex', 'ComplexTest', 'testDivideZeroZero', 234, 239]

    @Test
    public void testDivideZero() {
        Complex x = new Complex(3.0, 4.0);
        Complex z = x.divide(Complex.ZERO);
        // Assert.assertEquals(z, Complex.INF); // See MATH-657
        Assert.assertEquals(z, Complex.NaN);
    }

    @Test
    public void testDivideZeroZero() {
        Complex x = new Complex(0.0, 0.0);
        Complex z = x.divide(Complex.ZERO);
        Assert.assertEquals(z, Complex.NaN);
    }

[2140, 'src/test/java', 'org.apache.commons.math3.transform', 'FastCosineTransformerTest', 'testParameters', 217, 247]

[2140, 'src/test/java', 'org.apache.commons.math3.transform', 'FastSineTransformerTest', 'testParameters', 272, 299]

    @Test
    public void testParameters()
        throws Exception {
        UnivariateFunction f = new Sin();
        FastCosineTransformer transformer;
        transformer = new FastCosineTransformer(DctNormalization.STANDARD_DCT_I);

        try {
            // bad interval
            transformer.transform(f, 1, -1, 65, TransformType.FORWARD);
            Assert.fail("Expecting IllegalArgumentException - bad interval");
        } catch (IllegalArgumentException ex) {
            // expected
        }
        try {
            // bad samples number
            transformer.transform(f, -1, 1, 1, TransformType.FORWARD);
            Assert
                .fail("Expecting IllegalArgumentException - bad samples number");
        } catch (IllegalArgumentException ex) {
            // expected
        }
        try {
            // bad samples number
            transformer.transform(f, -1, 1, 64, TransformType.FORWARD);
            Assert
                .fail("Expecting IllegalArgumentException - bad samples number");
        } catch (IllegalArgumentException ex) {
            // expected
        }
    }

    @Test
    public void testParameters() throws Exception {
        UnivariateFunction f = new Sin();
        FastSineTransformer transformer;
        transformer = new FastSineTransformer(DstNormalization.STANDARD_DST_I);

        try {
            // bad interval
            transformer.transform(f, 1, -1, 64, TransformType.FORWARD);
            Assert.fail("Expecting IllegalArgumentException - bad interval");
        } catch (IllegalArgumentException ex) {
            // expected
        }
        try {
            // bad samples number
            transformer.transform(f, -1, 1, 0, TransformType.FORWARD);
            Assert.fail("Expecting IllegalArgumentException - bad samples number");
        } catch (IllegalArgumentException ex) {
            // expected
        }
        try {
            // bad samples number
            transformer.transform(f, -1, 1, 100, TransformType.FORWARD);
            Assert.fail("Expecting IllegalArgumentException - bad samples number");
        } catch (IllegalArgumentException ex) {
            // expected
        }
    }

[2148, 'src/test/java', 'org.apache.commons.math3.geometry.euclidean.threed', 'FieldRotationDSTest', 'testVectorTwoPairs', 405, 443]

[2148, 'src/test/java', 'org.apache.commons.math3.geometry.euclidean.threed', 'FieldRotationDfpTest', 'testVectorTwoPairs', 245, 283]

    @Test
    public void testVectorTwoPairs() throws MathArithmeticException {

        FieldVector3D<DerivativeStructure> u1 = createVector(3, 0, 0);
        FieldVector3D<DerivativeStructure> u2 = createVector(0, 5, 0);
        FieldVector3D<DerivativeStructure> v1 = createVector(0, 0, 2);
        FieldVector3D<DerivativeStructure> v2 = createVector(-2, 0, 2);
        FieldRotation<DerivativeStructure> r = new FieldRotation<DerivativeStructure>(u1, u2, v1, v2);
        checkVector(r.applyTo(createVector(1, 0, 0)), createVector(0, 0, 1));
        checkVector(r.applyTo(createVector(0, 1, 0)), createVector(-1, 0, 0));

        r = new FieldRotation<DerivativeStructure>(u1, u2, u1.negate(), u2.negate());
        FieldVector3D<DerivativeStructure> axis = r.getAxis(RotationConvention.VECTOR_OPERATOR);
        if (FieldVector3D.dotProduct(axis, createVector(0, 0, 1)).getReal() > 0) {
            checkVector(axis, createVector(0, 0, 1));
        } else {
            checkVector(axis, createVector(0, 0, -1));
        }
        checkAngle(r.getAngle(), FastMath.PI);

        double sqrt = FastMath.sqrt(2) / 2;
        r = new FieldRotation<DerivativeStructure>(createVector(1, 0, 0),  createVector(0, 1, 0),
                           createVector(0.5, 0.5,  sqrt),
                           createVector(0.5, 0.5, -sqrt));
        checkRotationDS(r, sqrt, 0.5, 0.5, 0);

        r = new FieldRotation<DerivativeStructure>(u1, u2, u1, FieldVector3D.crossProduct(u1, u2));
        checkRotationDS(r, sqrt, -sqrt, 0, 0);

        checkRotationDS(new FieldRotation<DerivativeStructure>(u1, u2, u1, u2), 1, 0, 0, 0);

        try {
            new FieldRotation<DerivativeStructure>(u1, u2, createVector(0, 0, 0), v2);
            Assert.fail("an exception should have been thrown");
        } catch (MathArithmeticException e) {
            // expected behavior
        }

    }

    @Test
    public void testVectorTwoPairs() throws MathArithmeticException {

        FieldVector3D<Dfp> u1 = createVector(3, 0, 0);
        FieldVector3D<Dfp> u2 = createVector(0, 5, 0);
        FieldVector3D<Dfp> v1 = createVector(0, 0, 2);
        FieldVector3D<Dfp> v2 = createVector(-2, 0, 2);
        FieldRotation<Dfp> r = new FieldRotation<Dfp>(u1, u2, v1, v2);
        checkVector(r.applyTo(createVector(1, 0, 0)), createVector(0, 0, 1));
        checkVector(r.applyTo(createVector(0, 1, 0)), createVector(-1, 0, 0));

        r = new FieldRotation<Dfp>(u1, u2, u1.negate(), u2.negate());
        FieldVector3D<Dfp> axis = r.getAxis(RotationConvention.VECTOR_OPERATOR);
        if (FieldVector3D.dotProduct(axis, createVector(0, 0, 1)).getReal() > 0) {
            checkVector(axis, createVector(0, 0, 1));
        } else {
            checkVector(axis, createVector(0, 0, -1));
        }
        checkAngle(r.getAngle(), FastMath.PI);

        double sqrt = FastMath.sqrt(2) / 2;
        r = new FieldRotation<Dfp>(createVector(1, 0, 0),  createVector(0, 1, 0),
                           createVector(0.5, 0.5,  sqrt),
                           createVector(0.5, 0.5, -sqrt));
        checkRotationDS(r, sqrt, 0.5, 0.5, 0);

        r = new FieldRotation<Dfp>(u1, u2, u1, FieldVector3D.crossProduct(u1, u2));
        checkRotationDS(r, sqrt, -sqrt, 0, 0);

        checkRotationDS(new FieldRotation<Dfp>(u1, u2, u1, u2), 1, 0, 0, 0);

        try {
            new FieldRotation<Dfp>(u1, u2, createVector(0, 0, 0), v2);
            Assert.fail("an exception should have been thrown");
        } catch (MathArithmeticException e) {
            // expected behavior
        }

    }

[2163, 'src/test/java', 'org.apache.commons.math3.linear', 'Array2DRowRealMatrixTest', 'testTrace', 280, 291]

[2163, 'src/test/java', 'org.apache.commons.math3.linear', 'BlockRealMatrixTest', 'testTrace', 305, 316]

    @Test
    public void testTrace() {
        RealMatrix m = new Array2DRowRealMatrix(id);
        Assert.assertEquals("identity trace",3d,m.getTrace(),entryTolerance);
        m = new Array2DRowRealMatrix(testData2);
        try {
            m.getTrace();
            Assert.fail("Expecting NonSquareMatrixException");
        } catch (NonSquareMatrixException ex) {
            // ignored
        }
    }

    @Test
    public void testTrace() {
        RealMatrix m = new BlockRealMatrix(id);
        Assert.assertEquals("identity trace",3d,m.getTrace(),entryTolerance);
        m = new BlockRealMatrix(testData2);
        try {
            m.getTrace();
            Assert.fail("Expecting NonSquareMatrixException");
        } catch (NonSquareMatrixException ex) {
            // ignored
        }
    }

[2173, 'src/test/java', 'org.apache.commons.math3.analysis.integration.gauss', 'LegendreHighPrecisionTest', 'testInverse', 44, 57]

[2173, 'src/test/java', 'org.apache.commons.math3.analysis.integration.gauss', 'LegendreTest', 'testInverse', 44, 57]

    @Test
    public void testInverse() {
        final UnivariateFunction inv = new Inverse();
        final UnivariateFunction log = new Log();

        final double lo = 12.34;
        final double hi = 456.78;

        final GaussIntegrator integrator = factory.legendreHighPrecision(60, lo, hi);
        final double s = integrator.integrate(inv);
        final double expected = log.value(hi) - log.value(lo);
        // System.out.println("s=" + s + " e=" + expected);
        Assert.assertEquals(expected, s, 1e-15);
    }

    @Test
    public void testInverse() {
        final UnivariateFunction inv = new Inverse();
        final UnivariateFunction log = new Log();

        final double lo = 12.34;
        final double hi = 456.78;

        final GaussIntegrator integrator = factory.legendre(60, lo, hi);
        final double s = integrator.integrate(inv);
        final double expected = log.value(hi) - log.value(lo);
        // System.out.println("s=" + s + " e=" + expected);
        Assert.assertEquals(expected, s, 1e-14);
    }

[2175, 'src/test/java', 'org.apache.commons.math3.linear', 'RealMatrixFormatAbstractTest', 'testSimpleWithDecimals', 54, 67]

[2175, 'src/test/java', 'org.apache.commons.math3.linear', 'RealMatrixFormatAbstractTest', 'testSimpleWithDecimalsTrunc', 69, 83]

    @Test
    public void testSimpleWithDecimals() {
        RealMatrix m = MatrixUtils.createRealMatrix(new double[][] {{1.23, 1.43, 1.63}, {2.46, 2.46, 2.66}});
        String expected =
            "{{1"    + getDecimalCharacter() +
            "23,1" + getDecimalCharacter() +
            "43,1" + getDecimalCharacter() +
            "63},{2" + getDecimalCharacter() +
            "46,2" + getDecimalCharacter() +
            "46,2" + getDecimalCharacter() +
            "66}}";
        String actual = realMatrixFormat.format(m);
        Assert.assertEquals(expected, actual);
    }

    @Test
    public void testSimpleWithDecimalsTrunc() {
        RealMatrix m = MatrixUtils.createRealMatrix(new double[][] {{1.232323232323, 1.43, 1.63},
                                                                    {2.46, 2.46, 2.666666666666}});
        String expected =
                "{{1"    + getDecimalCharacter() +
                "2323232323,1" + getDecimalCharacter() +
                "43,1" + getDecimalCharacter() +
                "63},{2" + getDecimalCharacter() +
                "46,2" + getDecimalCharacter() +
                "46,2" + getDecimalCharacter() +
                "6666666667}}";
        String actual = realMatrixFormat.format(m);
        Assert.assertEquals(expected, actual);
    }

[2181, 'src/test/java', 'org.apache.commons.math3.geometry.euclidean.threed', 'FieldRotationDSTest', 'testMatrix', 445, 589]

[2181, 'src/test/java', 'org.apache.commons.math3.geometry.euclidean.threed', 'FieldRotationDfpTest', 'testMatrix', 285, 429]

    @Test
    public void testMatrix()
            throws NotARotationMatrixException {

        try {
            createRotation(new double[][] {
                { 0.0, 1.0, 0.0 },
                { 1.0, 0.0, 0.0 }
            }, 1.0e-7);
            Assert.fail("Expecting NotARotationMatrixException");
        } catch (NotARotationMatrixException nrme) {
            // expected behavior
        }

        try {
            createRotation(new double[][] {
                {  0.445888,  0.797184, -0.407040 },
                {  0.821760, -0.184320,  0.539200 },
                { -0.354816,  0.574912,  0.737280 }
            }, 1.0e-7);
            Assert.fail("Expecting NotARotationMatrixException");
        } catch (NotARotationMatrixException nrme) {
            // expected behavior
        }

        try {
            createRotation(new double[][] {
                {  0.4,  0.8, -0.4 },
                { -0.4,  0.6,  0.7 },
                {  0.8, -0.2,  0.5 }
            }, 1.0e-15);
            Assert.fail("Expecting NotARotationMatrixException");
        } catch (NotARotationMatrixException nrme) {
            // expected behavior
        }

        checkRotationDS(createRotation(new double[][] {
            {  0.445888,  0.797184, -0.407040 },
            { -0.354816,  0.574912,  0.737280 },
            {  0.821760, -0.184320,  0.539200 }
        }, 1.0e-10),
        0.8, 0.288, 0.384, 0.36);

        checkRotationDS(createRotation(new double[][] {
            {  0.539200,  0.737280,  0.407040 },
            {  0.184320, -0.574912,  0.797184 },
            {  0.821760, -0.354816, -0.445888 }
        }, 1.0e-10),
        0.36, 0.8, 0.288, 0.384);

        checkRotationDS(createRotation(new double[][] {
            { -0.445888,  0.797184, -0.407040 },
            {  0.354816,  0.574912,  0.737280 },
            {  0.821760,  0.184320, -0.539200 }
        }, 1.0e-10),
        0.384, 0.36, 0.8, 0.288);

        checkRotationDS(createRotation(new double[][] {
            { -0.539200,  0.737280,  0.407040 },
            { -0.184320, -0.574912,  0.797184 },
            {  0.821760,  0.354816,  0.445888 }
        }, 1.0e-10),
        0.288, 0.384, 0.36, 0.8);

        double[][] m1 = { { 0.0, 1.0, 0.0 },
            { 0.0, 0.0, 1.0 },
            { 1.0, 0.0, 0.0 } };
        FieldRotation<DerivativeStructure> r = createRotation(m1, 1.0e-7);
        checkVector(r.applyTo(createVector(1, 0, 0)), createVector(0, 0, 1));
        checkVector(r.applyTo(createVector(0, 1, 0)), createVector(1, 0, 0));
        checkVector(r.applyTo(createVector(0, 0, 1)), createVector(0, 1, 0));

        double[][] m2 = { { 0.83203, -0.55012, -0.07139 },
            { 0.48293,  0.78164, -0.39474 },
            { 0.27296,  0.29396,  0.91602 } };
        r = createRotation(m2, 1.0e-12);

        DerivativeStructure[][] m3 = r.getMatrix();
        double d00 = m2[0][0] - m3[0][0].getReal();
        double d01 = m2[0][1] - m3[0][1].getReal();
        double d02 = m2[0][2] - m3[0][2].getReal();
        double d10 = m2[1][0] - m3[1][0].getReal();
        double d11 = m2[1][1] - m3[1][1].getReal();
        double d12 = m2[1][2] - m3[1][2].getReal();
        double d20 = m2[2][0] - m3[2][0].getReal();
        double d21 = m2[2][1] - m3[2][1].getReal();
        double d22 = m2[2][2] - m3[2][2].getReal();

        Assert.assertTrue(FastMath.abs(d00) < 6.0e-6);
        Assert.assertTrue(FastMath.abs(d01) < 6.0e-6);
        Assert.assertTrue(FastMath.abs(d02) < 6.0e-6);
        Assert.assertTrue(FastMath.abs(d10) < 6.0e-6);
        Assert.assertTrue(FastMath.abs(d11) < 6.0e-6);
        Assert.assertTrue(FastMath.abs(d12) < 6.0e-6);
        Assert.assertTrue(FastMath.abs(d20) < 6.0e-6);
        Assert.assertTrue(FastMath.abs(d21) < 6.0e-6);
        Assert.assertTrue(FastMath.abs(d22) < 6.0e-6);

        Assert.assertTrue(FastMath.abs(d00) > 4.0e-7);
        Assert.assertTrue(FastMath.abs(d01) > 4.0e-7);
        Assert.assertTrue(FastMath.abs(d02) > 4.0e-7);
        Assert.assertTrue(FastMath.abs(d10) > 4.0e-7);
        Assert.assertTrue(FastMath.abs(d11) > 4.0e-7);
        Assert.assertTrue(FastMath.abs(d12) > 4.0e-7);
        Assert.assertTrue(FastMath.abs(d20) > 4.0e-7);
        Assert.assertTrue(FastMath.abs(d21) > 4.0e-7);
        Assert.assertTrue(FastMath.abs(d22) > 4.0e-7);

        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                double m3tm3 = m3[i][0].getReal() * m3[j][0].getReal() +
                               m3[i][1].getReal() * m3[j][1].getReal() +
                               m3[i][2].getReal() * m3[j][2].getReal();
                if (i == j) {
                    Assert.assertTrue(FastMath.abs(m3tm3 - 1.0) < 1.0e-10);
                } else {
                    Assert.assertTrue(FastMath.abs(m3tm3) < 1.0e-10);
                }
            }
        }

        checkVector(r.applyTo(createVector(1, 0, 0)),
                    new FieldVector3D<DerivativeStructure>(m3[0][0], m3[1][0], m3[2][0]));
        checkVector(r.applyTo(createVector(0, 1, 0)),
                    new FieldVector3D<DerivativeStructure>(m3[0][1], m3[1][1], m3[2][1]));
        checkVector(r.applyTo(createVector(0, 0, 1)),
                    new FieldVector3D<DerivativeStructure>(m3[0][2], m3[1][2], m3[2][2]));

        double[][] m4 = { { 1.0,  0.0,  0.0 },
            { 0.0, -1.0,  0.0 },
            { 0.0,  0.0, -1.0 } };
        r = createRotation(m4, 1.0e-7);
        checkAngle(r.getAngle(), FastMath.PI);

        try {
            double[][] m5 = { { 0.0, 0.0, 1.0 },
                { 0.0, 1.0, 0.0 },
                { 1.0, 0.0, 0.0 } };
            r = createRotation(m5, 1.0e-7);
            Assert.fail("got " + r + ", should have caught an exception");
        } catch (NotARotationMatrixException e) {
            // expected
        }

    }

    @Test
    public void testMatrix()
            throws NotARotationMatrixException {

        try {
            createRotation(new double[][] {
                { 0.0, 1.0, 0.0 },
                { 1.0, 0.0, 0.0 }
            }, 1.0e-7);
            Assert.fail("Expecting NotARotationMatrixException");
        } catch (NotARotationMatrixException nrme) {
            // expected behavior
        }

        try {
            createRotation(new double[][] {
                {  0.445888,  0.797184, -0.407040 },
                {  0.821760, -0.184320,  0.539200 },
                { -0.354816,  0.574912,  0.737280 }
            }, 1.0e-7);
            Assert.fail("Expecting NotARotationMatrixException");
        } catch (NotARotationMatrixException nrme) {
            // expected behavior
        }

        try {
            createRotation(new double[][] {
                {  0.4,  0.8, -0.4 },
                { -0.4,  0.6,  0.7 },
                {  0.8, -0.2,  0.5 }
            }, 1.0e-15);
            Assert.fail("Expecting NotARotationMatrixException");
        } catch (NotARotationMatrixException nrme) {
            // expected behavior
        }

        checkRotationDS(createRotation(new double[][] {
            {  0.445888,  0.797184, -0.407040 },
            { -0.354816,  0.574912,  0.737280 },
            {  0.821760, -0.184320,  0.539200 }
        }, 1.0e-10),
        0.8, 0.288, 0.384, 0.36);

        checkRotationDS(createRotation(new double[][] {
            {  0.539200,  0.737280,  0.407040 },
            {  0.184320, -0.574912,  0.797184 },
            {  0.821760, -0.354816, -0.445888 }
        }, 1.0e-10),
        0.36, 0.8, 0.288, 0.384);

        checkRotationDS(createRotation(new double[][] {
            { -0.445888,  0.797184, -0.407040 },
            {  0.354816,  0.574912,  0.737280 },
            {  0.821760,  0.184320, -0.539200 }
        }, 1.0e-10),
        0.384, 0.36, 0.8, 0.288);

        checkRotationDS(createRotation(new double[][] {
            { -0.539200,  0.737280,  0.407040 },
            { -0.184320, -0.574912,  0.797184 },
            {  0.821760,  0.354816,  0.445888 }
        }, 1.0e-10),
        0.288, 0.384, 0.36, 0.8);

        double[][] m1 = { { 0.0, 1.0, 0.0 },
            { 0.0, 0.0, 1.0 },
            { 1.0, 0.0, 0.0 } };
        FieldRotation<Dfp> r = createRotation(m1, 1.0e-7);
        checkVector(r.applyTo(createVector(1, 0, 0)), createVector(0, 0, 1));
        checkVector(r.applyTo(createVector(0, 1, 0)), createVector(1, 0, 0));
        checkVector(r.applyTo(createVector(0, 0, 1)), createVector(0, 1, 0));

        double[][] m2 = { { 0.83203, -0.55012, -0.07139 },
            { 0.48293,  0.78164, -0.39474 },
            { 0.27296,  0.29396,  0.91602 } };
        r = createRotation(m2, 1.0e-12);

        Dfp[][] m3 = r.getMatrix();
        double d00 = m2[0][0] - m3[0][0].getReal();
        double d01 = m2[0][1] - m3[0][1].getReal();
        double d02 = m2[0][2] - m3[0][2].getReal();
        double d10 = m2[1][0] - m3[1][0].getReal();
        double d11 = m2[1][1] - m3[1][1].getReal();
        double d12 = m2[1][2] - m3[1][2].getReal();
        double d20 = m2[2][0] - m3[2][0].getReal();
        double d21 = m2[2][1] - m3[2][1].getReal();
        double d22 = m2[2][2] - m3[2][2].getReal();

        Assert.assertTrue(FastMath.abs(d00) < 6.0e-6);
        Assert.assertTrue(FastMath.abs(d01) < 6.0e-6);
        Assert.assertTrue(FastMath.abs(d02) < 6.0e-6);
        Assert.assertTrue(FastMath.abs(d10) < 6.0e-6);
        Assert.assertTrue(FastMath.abs(d11) < 6.0e-6);
        Assert.assertTrue(FastMath.abs(d12) < 6.0e-6);
        Assert.assertTrue(FastMath.abs(d20) < 6.0e-6);
        Assert.assertTrue(FastMath.abs(d21) < 6.0e-6);
        Assert.assertTrue(FastMath.abs(d22) < 6.0e-6);

        Assert.assertTrue(FastMath.abs(d00) > 4.0e-7);
        Assert.assertTrue(FastMath.abs(d01) > 4.0e-7);
        Assert.assertTrue(FastMath.abs(d02) > 4.0e-7);
        Assert.assertTrue(FastMath.abs(d10) > 4.0e-7);
        Assert.assertTrue(FastMath.abs(d11) > 4.0e-7);
        Assert.assertTrue(FastMath.abs(d12) > 4.0e-7);
        Assert.assertTrue(FastMath.abs(d20) > 4.0e-7);
        Assert.assertTrue(FastMath.abs(d21) > 4.0e-7);
        Assert.assertTrue(FastMath.abs(d22) > 4.0e-7);

        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                double m3tm3 = m3[i][0].getReal() * m3[j][0].getReal() +
                               m3[i][1].getReal() * m3[j][1].getReal() +
                               m3[i][2].getReal() * m3[j][2].getReal();
                if (i == j) {
                    Assert.assertTrue(FastMath.abs(m3tm3 - 1.0) < 1.0e-10);
                } else {
                    Assert.assertTrue(FastMath.abs(m3tm3) < 1.0e-10);
                }
            }
        }

        checkVector(r.applyTo(createVector(1, 0, 0)),
                    new FieldVector3D<Dfp>(m3[0][0], m3[1][0], m3[2][0]));
        checkVector(r.applyTo(createVector(0, 1, 0)),
                    new FieldVector3D<Dfp>(m3[0][1], m3[1][1], m3[2][1]));
        checkVector(r.applyTo(createVector(0, 0, 1)),
                    new FieldVector3D<Dfp>(m3[0][2], m3[1][2], m3[2][2]));

        double[][] m4 = { { 1.0,  0.0,  0.0 },
            { 0.0, -1.0,  0.0 },
            { 0.0,  0.0, -1.0 } };
        r = createRotation(m4, 1.0e-7);
        checkAngle(r.getAngle(), FastMath.PI);

        try {
            double[][] m5 = { { 0.0, 0.0, 1.0 },
                { 0.0, 1.0, 0.0 },
                { 1.0, 0.0, 0.0 } };
            r = createRotation(m5, 1.0e-7);
            Assert.fail("got " + r + ", should have caught an exception");
        } catch (NotARotationMatrixException e) {
            // expected
        }

    }

[2191, 'src/test/java', 'org.apache.commons.math3.distribution', 'ChiSquaredDistributionTest', 'testMoments', 123, 135]

[2191, 'src/test/java', 'org.apache.commons.math3.distribution', 'PoissonDistributionTest', 'testMoments', 233, 245]

    @Test
    public void testMoments() {
        final double tol = 1e-9;
        ChiSquaredDistribution dist;

        dist = new ChiSquaredDistribution(1500);
        Assert.assertEquals(dist.getNumericalMean(), 1500, tol);
        Assert.assertEquals(dist.getNumericalVariance(), 3000, tol);

        dist = new ChiSquaredDistribution(1.12);
        Assert.assertEquals(dist.getNumericalMean(), 1.12, tol);
        Assert.assertEquals(dist.getNumericalVariance(), 2.24, tol);
    }

    @Test
    public void testMoments() {
        final double tol = 1e-9;
        PoissonDistribution dist;

        dist = new PoissonDistribution(1);
        Assert.assertEquals(dist.getNumericalMean(), 1, tol);
        Assert.assertEquals(dist.getNumericalVariance(), 1, tol);

        dist = new PoissonDistribution(11.23);
        Assert.assertEquals(dist.getNumericalMean(), 11.23, tol);
        Assert.assertEquals(dist.getNumericalVariance(), 11.23, tol);
    }

[2192, 'src/test/java', 'org.apache.commons.math3.stat.correlation', 'KendallsCorrelationTest', 'testSimpleOrdered', 81, 91]

[2192, 'src/test/java', 'org.apache.commons.math3.stat.correlation', 'KendallsCorrelationTest', 'testSimpleOrderedPowerOf2', 105, 115]

    @Test
    public void testSimpleOrdered() {
        final int length = 10;
        final double[] xArray = new double[length];
        final double[] yArray = new double[length];
        for (int i = 0; i < length; i++) {
            xArray[i] = i;
            yArray[i] = i;
        }
        Assert.assertEquals(1.0, correlation.correlation(xArray, yArray), Double.MIN_VALUE);
    }

    @Test
    public void testSimpleOrderedPowerOf2() {
        final int length = 16;
        final double[] xArray = new double[length];
        final double[] yArray = new double[length];
        for (int i = 0; i < length; i++) {
            xArray[i] = i;
            yArray[i] = i;
        }
        Assert.assertEquals(1.0, correlation.correlation(xArray, yArray), Double.MIN_VALUE);
    }

[2209, 'src/test/java', 'org.apache.commons.math3.geometry.euclidean.threed', 'FieldRotationDSTest', 'testIdentity', 35, 56]

[2209, 'src/test/java', 'org.apache.commons.math3.geometry.euclidean.threed', 'FieldRotationDfpTest', 'testIdentity', 34, 55]

    @Test
    public void testIdentity() {

        FieldRotation<DerivativeStructure> r = createRotation(1, 0, 0, 0, false);
        checkVector(r.applyTo(createVector(1, 0, 0)), createVector(1, 0, 0));
        checkVector(r.applyTo(createVector(0, 1, 0)), createVector(0, 1, 0));
        checkVector(r.applyTo(createVector(0, 0, 1)), createVector(0, 0, 1));
        checkAngle(r.getAngle(), 0);

        r = createRotation(-1, 0, 0, 0, false);
        checkVector(r.applyTo(createVector(1, 0, 0)), createVector(1, 0, 0));
        checkVector(r.applyTo(createVector(0, 1, 0)), createVector(0, 1, 0));
        checkVector(r.applyTo(createVector(0, 0, 1)), createVector(0, 0, 1));
        checkAngle(r.getAngle(), 0);

        r = createRotation(42, 0, 0, 0, true);
        checkVector(r.applyTo(createVector(1, 0, 0)), createVector(1, 0, 0));
        checkVector(r.applyTo(createVector(0, 1, 0)), createVector(0, 1, 0));
        checkVector(r.applyTo(createVector(0, 0, 1)), createVector(0, 0, 1));
        checkAngle(r.getAngle(), 0);

    }

    @Test
    public void testIdentity() {

        FieldRotation<Dfp> r = createRotation(1, 0, 0, 0, false);
        checkVector(r.applyTo(createVector(1, 0, 0)), createVector(1, 0, 0));
        checkVector(r.applyTo(createVector(0, 1, 0)), createVector(0, 1, 0));
        checkVector(r.applyTo(createVector(0, 0, 1)), createVector(0, 0, 1));
        checkAngle(r.getAngle(), 0);

        r = createRotation(-1, 0, 0, 0, false);
        checkVector(r.applyTo(createVector(1, 0, 0)), createVector(1, 0, 0));
        checkVector(r.applyTo(createVector(0, 1, 0)), createVector(0, 1, 0));
        checkVector(r.applyTo(createVector(0, 0, 1)), createVector(0, 0, 1));
        checkAngle(r.getAngle(), 0);

        r = createRotation(42, 0, 0, 0, true);
        checkVector(r.applyTo(createVector(1, 0, 0)), createVector(1, 0, 0));
        checkVector(r.applyTo(createVector(0, 1, 0)), createVector(0, 1, 0));
        checkVector(r.applyTo(createVector(0, 0, 1)), createVector(0, 0, 1));
        checkAngle(r.getAngle(), 0);

    }

[2213, 'src/test/java', 'org.apache.commons.math3.optim.nonlinear.scalar', 'MultivariateFunctionMappingAdapterTest', 'testStartSimplexInsideRange', 30, 55]

[2213, 'src/test/java', 'org.apache.commons.math3.optim.nonlinear.scalar', 'MultivariateFunctionMappingAdapterTest', 'testOptimumOutsideRange', 57, 82]

    @Test
    public void testStartSimplexInsideRange() {
        final BiQuadratic biQuadratic = new BiQuadratic(2.0, 2.5, 1.0, 3.0, 2.0, 3.0);
        final MultivariateFunctionMappingAdapter wrapped
            = new MultivariateFunctionMappingAdapter(biQuadratic,
                                                     biQuadratic.getLower(),
                                                     biQuadratic.getUpper());

        SimplexOptimizer optimizer = new SimplexOptimizer(1e-10, 1e-30);
        final AbstractSimplex simplex = new NelderMeadSimplex(new double[][] {
                wrapped.boundedToUnbounded(new double[] { 1.5, 2.75 }),
                wrapped.boundedToUnbounded(new double[] { 1.5, 2.95 }),
                wrapped.boundedToUnbounded(new double[] { 1.7, 2.90 })
            });

        final PointValuePair optimum
            = optimizer.optimize(new MaxEval(300),
                                 new ObjectiveFunction(wrapped),
                                 simplex,
                                 GoalType.MINIMIZE,
                                 new InitialGuess(wrapped.boundedToUnbounded(new double[] { 1.5, 2.25 })));
        final double[] bounded = wrapped.unboundedToBounded(optimum.getPoint());

        Assert.assertEquals(biQuadratic.getBoundedXOptimum(), bounded[0], 2e-7);
        Assert.assertEquals(biQuadratic.getBoundedYOptimum(), bounded[1], 2e-7);
    }

    @Test
    public void testOptimumOutsideRange() {
        final BiQuadratic biQuadratic = new BiQuadratic(4.0, 0.0, 1.0, 3.0, 2.0, 3.0);
        final MultivariateFunctionMappingAdapter wrapped
            = new MultivariateFunctionMappingAdapter(biQuadratic,
                                                     biQuadratic.getLower(),
                                                     biQuadratic.getUpper());

        SimplexOptimizer optimizer = new SimplexOptimizer(1e-10, 1e-30);
        final AbstractSimplex simplex = new NelderMeadSimplex(new double[][] {
                wrapped.boundedToUnbounded(new double[] { 1.5, 2.75 }),
                wrapped.boundedToUnbounded(new double[] { 1.5, 2.95 }),
                wrapped.boundedToUnbounded(new double[] { 1.7, 2.90 })
            });

        final PointValuePair optimum
            = optimizer.optimize(new MaxEval(100),
                                 new ObjectiveFunction(wrapped),
                                 simplex,
                                 GoalType.MINIMIZE,
                                 new InitialGuess(wrapped.boundedToUnbounded(new double[] { 1.5, 2.25 })));
        final double[] bounded = wrapped.unboundedToBounded(optimum.getPoint());

        Assert.assertEquals(biQuadratic.getBoundedXOptimum(), bounded[0], 2e-7);
        Assert.assertEquals(biQuadratic.getBoundedYOptimum(), bounded[1], 2e-7);
    }

[2222, 'src/test/java', 'org.apache.commons.math3.linear', 'RealVectorAbstractTest', 'testMap', 884, 890]

[2222, 'src/test/java', 'org.apache.commons.math3.linear', 'RealVectorAbstractTest', 'testMapToSelf', 892, 898]

    @Test
    public void testMap() {
        final UnivariateFunction[] functions = createFunctions();
        for (UnivariateFunction f : functions) {
            doTestMapFunction(f, false);
        }
    }

    @Test
    public void testMapToSelf() {
        final UnivariateFunction[] functions = createFunctions();
        for (UnivariateFunction f : functions) {
            doTestMapFunction(f, true);
        }
    }

[2226, 'src/test/java', 'org.apache.commons.math3.stat.inference', 'TTestTest', 'testOneSampleTTest', 123, 159]

[2226, 'src/test/java', 'org.apache.commons.math3.stat.inference', 'TestUtilsTest', 'testOneSampleTTest', 276, 312]

    @Test
    public void testOneSampleTTest() {
        double[] oneSidedP =
            {2d, 0d, 6d, 6d, 3d, 3d, 2d, 3d, -6d, 6d, 6d, 6d, 3d, 0d, 1d, 1d, 0d, 2d, 3d, 3d };
        SummaryStatistics oneSidedPStats = new SummaryStatistics();
        for (int i = 0; i < oneSidedP.length; i++) {
            oneSidedPStats.addValue(oneSidedP[i]);
        }
        // Target comparison values computed using R version 1.8.1 (Linux version)
        Assert.assertEquals("one sample t stat", 3.86485535541,
                testStatistic.t(0d, oneSidedP), 10E-10);
        Assert.assertEquals("one sample t stat", 3.86485535541,
                testStatistic.t(0d, oneSidedPStats),1E-10);
        Assert.assertEquals("one sample p value", 0.000521637019637,
                testStatistic.tTest(0d, oneSidedP) / 2d, 10E-10);
        Assert.assertEquals("one sample p value", 0.000521637019637,
                testStatistic.tTest(0d, oneSidedPStats) / 2d, 10E-5);
        Assert.assertTrue("one sample t-test reject", testStatistic.tTest(0d, oneSidedP, 0.01));
        Assert.assertTrue("one sample t-test reject", testStatistic.tTest(0d, oneSidedPStats, 0.01));
        Assert.assertTrue("one sample t-test accept", !testStatistic.tTest(0d, oneSidedP, 0.0001));
        Assert.assertTrue("one sample t-test accept", !testStatistic.tTest(0d, oneSidedPStats, 0.0001));

        try {
            testStatistic.tTest(0d, oneSidedP, 95);
            Assert.fail("alpha out of range, OutOfRangeException expected");
        } catch (OutOfRangeException ex) {
            // expected
        }

        try {
            testStatistic.tTest(0d, oneSidedPStats, 95);
            Assert.fail("alpha out of range, OutOfRangeException expected");
        } catch (OutOfRangeException ex) {
            // expected
        }

    }

    @Test
    public void testOneSampleTTest() {
        double[] oneSidedP =
            {2d, 0d, 6d, 6d, 3d, 3d, 2d, 3d, -6d, 6d, 6d, 6d, 3d, 0d, 1d, 1d, 0d, 2d, 3d, 3d };
        SummaryStatistics oneSidedPStats = new SummaryStatistics();
        for (int i = 0; i < oneSidedP.length; i++) {
            oneSidedPStats.addValue(oneSidedP[i]);
        }
        // Target comparison values computed using R version 1.8.1 (Linux version)
        Assert.assertEquals("one sample t stat", 3.86485535541,
                TestUtils.t(0d, oneSidedP), 10E-10);
        Assert.assertEquals("one sample t stat", 3.86485535541,
                TestUtils.t(0d, oneSidedPStats),1E-10);
        Assert.assertEquals("one sample p value", 0.000521637019637,
                TestUtils.tTest(0d, oneSidedP) / 2d, 10E-10);
        Assert.assertEquals("one sample p value", 0.000521637019637,
                TestUtils.tTest(0d, oneSidedPStats) / 2d, 10E-5);
        Assert.assertTrue("one sample t-test reject", TestUtils.tTest(0d, oneSidedP, 0.01));
        Assert.assertTrue("one sample t-test reject", TestUtils.tTest(0d, oneSidedPStats, 0.01));
        Assert.assertTrue("one sample t-test accept", !TestUtils.tTest(0d, oneSidedP, 0.0001));
        Assert.assertTrue("one sample t-test accept", !TestUtils.tTest(0d, oneSidedPStats, 0.0001));

        try {
            TestUtils.tTest(0d, oneSidedP, 95);
            Assert.fail("alpha out of range, OutOfRangeException expected");
        } catch (OutOfRangeException ex) {
            // expected
        }

        try {
            TestUtils.tTest(0d, oneSidedPStats, 95);
            Assert.fail("alpha out of range, OutOfRangeException expected");
        } catch (OutOfRangeException ex) {
            // expected
        }

    }

[2233, 'src/test/java', 'org.apache.commons.math3.ode.nonstiff', 'DormandPrince54IntegratorTest', 'testKepler', 233, 255]

[2233, 'src/test/java', 'org.apache.commons.math3.ode.nonstiff', 'DormandPrince853IntegratorTest', 'testKepler', 277, 299]

  @Test
  public void testKepler()
      throws DimensionMismatchException, NumberIsTooSmallException,
             MaxCountExceededException, NoBracketingException {

    final TestProblem3 pb  = new TestProblem3(0.9);
    double minStep = 0;
    double maxStep = pb.getFinalTime() - pb.getInitialTime();
    double scalAbsoluteTolerance = 1.0e-8;
    double scalRelativeTolerance = scalAbsoluteTolerance;

    FirstOrderIntegrator integ = new DormandPrince54Integrator(minStep, maxStep,
                                                               scalAbsoluteTolerance,
                                                               scalRelativeTolerance);
    integ.addStepHandler(new KeplerHandler(pb));
    integ.integrate(pb,
                    pb.getInitialTime(), pb.getInitialState(),
                    pb.getFinalTime(), new double[pb.getDimension()]);

    Assert.assertEquals(integ.getEvaluations(), pb.getCalls());
    Assert.assertTrue(pb.getCalls() < 2800);

  }

  @Test
  public void testKepler()
      throws DimensionMismatchException, NumberIsTooSmallException,
             MaxCountExceededException, NoBracketingException {

    final TestProblem3 pb  = new TestProblem3(0.9);
    double minStep = 0;
    double maxStep = pb.getFinalTime() - pb.getInitialTime();
    double scalAbsoluteTolerance = 1.0e-8;
    double scalRelativeTolerance = scalAbsoluteTolerance;

    FirstOrderIntegrator integ = new DormandPrince853Integrator(minStep, maxStep,
                                                                scalAbsoluteTolerance,
                                                                scalRelativeTolerance);
    integ.addStepHandler(new KeplerHandler(pb));
    integ.integrate(pb,
                    pb.getInitialTime(), pb.getInitialState(),
                    pb.getFinalTime(), new double[pb.getDimension()]);

    Assert.assertEquals(integ.getEvaluations(), pb.getCalls());
    Assert.assertTrue(pb.getCalls() < 3300);

  }

[2302, 'src/test/java', 'org.apache.commons.math3.linear', 'BlockFieldMatrixTest', 'testWalk', 1242, 1327]

[2302, 'src/test/java', 'org.apache.commons.math3.linear', 'FieldMatrixImplTest', 'testWalk', 969, 1054]

    @Test
    public void testWalk() {
        int rows    = 150;
        int columns = 75;

        FieldMatrix<Fraction> m = new BlockFieldMatrix<Fraction>(FractionField.getInstance(), rows, columns);
        m.walkInRowOrder(new SetVisitor());
        GetVisitor getVisitor = new GetVisitor();
        m.walkInOptimizedOrder(getVisitor);
        Assert.assertEquals(rows * columns, getVisitor.getCount());

        m = new BlockFieldMatrix<Fraction>(FractionField.getInstance(), rows, columns);
        m.walkInRowOrder(new SetVisitor(), 1, rows - 2, 1, columns - 2);
        getVisitor = new GetVisitor();
        m.walkInOptimizedOrder(getVisitor, 1, rows - 2, 1, columns - 2);
        Assert.assertEquals((rows - 2) * (columns - 2), getVisitor.getCount());
        for (int i = 0; i < rows; ++i) {
            Assert.assertEquals(new Fraction(0), m.getEntry(i, 0));
            Assert.assertEquals(new Fraction(0), m.getEntry(i, columns - 1));
        }
        for (int j = 0; j < columns; ++j) {
            Assert.assertEquals(new Fraction(0), m.getEntry(0, j));
            Assert.assertEquals(new Fraction(0), m.getEntry(rows - 1, j));
        }

        m = new BlockFieldMatrix<Fraction>(FractionField.getInstance(), rows, columns);
        m.walkInColumnOrder(new SetVisitor());
        getVisitor = new GetVisitor();
        m.walkInOptimizedOrder(getVisitor);
        Assert.assertEquals(rows * columns, getVisitor.getCount());

        m = new BlockFieldMatrix<Fraction>(FractionField.getInstance(), rows, columns);
        m.walkInColumnOrder(new SetVisitor(), 1, rows - 2, 1, columns - 2);
        getVisitor = new GetVisitor();
        m.walkInOptimizedOrder(getVisitor, 1, rows - 2, 1, columns - 2);
        Assert.assertEquals((rows - 2) * (columns - 2), getVisitor.getCount());
        for (int i = 0; i < rows; ++i) {
            Assert.assertEquals(new Fraction(0), m.getEntry(i, 0));
            Assert.assertEquals(new Fraction(0), m.getEntry(i, columns - 1));
        }
        for (int j = 0; j < columns; ++j) {
            Assert.assertEquals(new Fraction(0), m.getEntry(0, j));
            Assert.assertEquals(new Fraction(0), m.getEntry(rows - 1, j));
        }

        m = new BlockFieldMatrix<Fraction>(FractionField.getInstance(), rows, columns);
        m.walkInOptimizedOrder(new SetVisitor());
        getVisitor = new GetVisitor();
        m.walkInRowOrder(getVisitor);
        Assert.assertEquals(rows * columns, getVisitor.getCount());

        m = new BlockFieldMatrix<Fraction>(FractionField.getInstance(), rows, columns);
        m.walkInOptimizedOrder(new SetVisitor(), 1, rows - 2, 1, columns - 2);
        getVisitor = new GetVisitor();
        m.walkInRowOrder(getVisitor, 1, rows - 2, 1, columns - 2);
        Assert.assertEquals((rows - 2) * (columns - 2), getVisitor.getCount());
        for (int i = 0; i < rows; ++i) {
            Assert.assertEquals(new Fraction(0), m.getEntry(i, 0));
            Assert.assertEquals(new Fraction(0), m.getEntry(i, columns - 1));
        }
        for (int j = 0; j < columns; ++j) {
            Assert.assertEquals(new Fraction(0), m.getEntry(0, j));
            Assert.assertEquals(new Fraction(0), m.getEntry(rows - 1, j));
        }

        m = new BlockFieldMatrix<Fraction>(FractionField.getInstance(), rows, columns);
        m.walkInOptimizedOrder(new SetVisitor());
        getVisitor = new GetVisitor();
        m.walkInColumnOrder(getVisitor);
        Assert.assertEquals(rows * columns, getVisitor.getCount());

        m = new BlockFieldMatrix<Fraction>(FractionField.getInstance(), rows, columns);
        m.walkInOptimizedOrder(new SetVisitor(), 1, rows - 2, 1, columns - 2);
        getVisitor = new GetVisitor();
        m.walkInColumnOrder(getVisitor, 1, rows - 2, 1, columns - 2);
        Assert.assertEquals((rows - 2) * (columns - 2), getVisitor.getCount());
        for (int i = 0; i < rows; ++i) {
            Assert.assertEquals(new Fraction(0), m.getEntry(i, 0));
            Assert.assertEquals(new Fraction(0), m.getEntry(i, columns - 1));
        }
        for (int j = 0; j < columns; ++j) {
            Assert.assertEquals(new Fraction(0), m.getEntry(0, j));
            Assert.assertEquals(new Fraction(0), m.getEntry(rows - 1, j));
        }

    }

    @Test
    public void testWalk() {
        int rows    = 150;
        int columns = 75;

        FieldMatrix<Fraction> m =
            new Array2DRowFieldMatrix<Fraction>(FractionField.getInstance(), rows, columns);
        m.walkInRowOrder(new SetVisitor());
        GetVisitor getVisitor = new GetVisitor();
        m.walkInOptimizedOrder(getVisitor);
        Assert.assertEquals(rows * columns, getVisitor.getCount());

        m = new Array2DRowFieldMatrix<Fraction>(FractionField.getInstance(), rows, columns);
        m.walkInRowOrder(new SetVisitor(), 1, rows - 2, 1, columns - 2);
        getVisitor = new GetVisitor();
        m.walkInOptimizedOrder(getVisitor, 1, rows - 2, 1, columns - 2);
        Assert.assertEquals((rows - 2) * (columns - 2), getVisitor.getCount());
        for (int i = 0; i < rows; ++i) {
            Assert.assertEquals(new Fraction(0), m.getEntry(i, 0));
            Assert.assertEquals(new Fraction(0), m.getEntry(i, columns - 1));
        }
        for (int j = 0; j < columns; ++j) {
            Assert.assertEquals(new Fraction(0), m.getEntry(0, j));
            Assert.assertEquals(new Fraction(0), m.getEntry(rows - 1, j));
        }

        m = new Array2DRowFieldMatrix<Fraction>(FractionField.getInstance(), rows, columns);
        m.walkInColumnOrder(new SetVisitor());
        getVisitor = new GetVisitor();
        m.walkInOptimizedOrder(getVisitor);
        Assert.assertEquals(rows * columns, getVisitor.getCount());

        m = new Array2DRowFieldMatrix<Fraction>(FractionField.getInstance(), rows, columns);
        m.walkInColumnOrder(new SetVisitor(), 1, rows - 2, 1, columns - 2);
        getVisitor = new GetVisitor();
        m.walkInOptimizedOrder(getVisitor, 1, rows - 2, 1, columns - 2);
        Assert.assertEquals((rows - 2) * (columns - 2), getVisitor.getCount());
        for (int i = 0; i < rows; ++i) {
            Assert.assertEquals(new Fraction(0), m.getEntry(i, 0));
            Assert.assertEquals(new Fraction(0), m.getEntry(i, columns - 1));
        }
        for (int j = 0; j < columns; ++j) {
            Assert.assertEquals(new Fraction(0), m.getEntry(0, j));
            Assert.assertEquals(new Fraction(0), m.getEntry(rows - 1, j));
        }

        m = new Array2DRowFieldMatrix<Fraction>(FractionField.getInstance(), rows, columns);
        m.walkInOptimizedOrder(new SetVisitor());
        getVisitor = new GetVisitor();
        m.walkInRowOrder(getVisitor);
        Assert.assertEquals(rows * columns, getVisitor.getCount());

        m = new Array2DRowFieldMatrix<Fraction>(FractionField.getInstance(), rows, columns);
        m.walkInOptimizedOrder(new SetVisitor(), 1, rows - 2, 1, columns - 2);
        getVisitor = new GetVisitor();
        m.walkInRowOrder(getVisitor, 1, rows - 2, 1, columns - 2);
        Assert.assertEquals((rows - 2) * (columns - 2), getVisitor.getCount());
        for (int i = 0; i < rows; ++i) {
            Assert.assertEquals(new Fraction(0), m.getEntry(i, 0));
            Assert.assertEquals(new Fraction(0), m.getEntry(i, columns - 1));
        }
        for (int j = 0; j < columns; ++j) {
            Assert.assertEquals(new Fraction(0), m.getEntry(0, j));
            Assert.assertEquals(new Fraction(0), m.getEntry(rows - 1, j));
        }

        m = new Array2DRowFieldMatrix<Fraction>(FractionField.getInstance(), rows, columns);
        m.walkInOptimizedOrder(new SetVisitor());
        getVisitor = new GetVisitor();
        m.walkInColumnOrder(getVisitor);
        Assert.assertEquals(rows * columns, getVisitor.getCount());

        m = new Array2DRowFieldMatrix<Fraction>(FractionField.getInstance(), rows, columns);
        m.walkInOptimizedOrder(new SetVisitor(), 1, rows - 2, 1, columns - 2);
        getVisitor = new GetVisitor();
        m.walkInColumnOrder(getVisitor, 1, rows - 2, 1, columns - 2);
        Assert.assertEquals((rows - 2) * (columns - 2), getVisitor.getCount());
        for (int i = 0; i < rows; ++i) {
            Assert.assertEquals(new Fraction(0), m.getEntry(i, 0));
            Assert.assertEquals(new Fraction(0), m.getEntry(i, columns - 1));
        }
        for (int j = 0; j < columns; ++j) {
            Assert.assertEquals(new Fraction(0), m.getEntry(0, j));
            Assert.assertEquals(new Fraction(0), m.getEntry(rows - 1, j));
        }
    }

[2323, 'src/test/java', 'org.apache.commons.math3.linear', 'DiagonalMatrixTest', 'testOperate', 235, 246]

[2323, 'src/test/java', 'org.apache.commons.math3.linear', 'DiagonalMatrixTest', 'testPreMultiply', 248, 259]

    @Test
    public void testOperate() {
        final double[] data = { -1.2, 3.4, 5 };
        final DiagonalMatrix diag = new DiagonalMatrix(data);
        final RealMatrix dense = new Array2DRowRealMatrix(diag.getData());

        final double[] v = { 6.7, 890.1, 23.4 };
        final double[] diagResult = diag.operate(v);
        final double[] denseResult = dense.operate(v);

        TestUtils.assertEquals(diagResult, denseResult, 0d);
    }

    @Test
    public void testPreMultiply() {
        final double[] data = { -1.2, 3.4, 5 };
        final DiagonalMatrix diag = new DiagonalMatrix(data);
        final RealMatrix dense = new Array2DRowRealMatrix(diag.getData());

        final double[] v = { 6.7, 890.1, 23.4 };
        final double[] diagResult = diag.preMultiply(v);
        final double[] denseResult = dense.preMultiply(v);

        TestUtils.assertEquals(diagResult, denseResult, 0d);
    }

[2333, 'src/test/java', 'org.apache.commons.math3.ode.nonstiff', 'EulerIntegratorTest', 'testDecreasingSteps', 55, 101]

[2333, 'src/test/java', 'org.apache.commons.math3.ode.nonstiff', 'MidpointIntegratorTest', 'testDecreasingSteps', 55, 101]

  @Test
  public void testDecreasingSteps()
      throws DimensionMismatchException, NumberIsTooSmallException,
             MaxCountExceededException, NoBracketingException {

      for (TestProblemAbstract pb : new TestProblemAbstract[] {
          new TestProblem1(), new TestProblem2(), new TestProblem3(),
          new TestProblem4(), new TestProblem5(), new TestProblem6()
      }) {

      double previousValueError = Double.NaN;
      double previousTimeError = Double.NaN;
      for (int i = 4; i < 8; ++i) {

        double step = (pb.getFinalTime() - pb.getInitialTime()) * FastMath.pow(2.0, -i);

        FirstOrderIntegrator integ = new EulerIntegrator(step);
        TestProblemHandler handler = new TestProblemHandler(pb, integ);
        integ.addStepHandler(handler);
        EventHandler[] functions = pb.getEventsHandlers();
        for (int l = 0; l < functions.length; ++l) {
          integ.addEventHandler(functions[l],
                                     Double.POSITIVE_INFINITY, 1.0e-6 * step, 1000);
        }
        double stopTime = integ.integrate(pb, pb.getInitialTime(), pb.getInitialState(),
                                          pb.getFinalTime(), new double[pb.getDimension()]);
        if (functions.length == 0) {
            Assert.assertEquals(pb.getFinalTime(), stopTime, 1.0e-10);
        }

        double valueError = handler.getMaximalValueError();
        if (i > 4) {
          Assert.assertTrue(valueError < FastMath.abs(previousValueError));
        }
        previousValueError = valueError;

        double timeError = handler.getMaximalTimeError();
        if (i > 4) {
          Assert.assertTrue(timeError <= FastMath.abs(previousTimeError));
        }
        previousTimeError = timeError;

      }

    }

  }

  @Test
  public void testDecreasingSteps()
      throws DimensionMismatchException, NumberIsTooSmallException,
             MaxCountExceededException, NoBracketingException {

      for (TestProblemAbstract pb : new TestProblemAbstract[] {
          new TestProblem1(), new TestProblem2(), new TestProblem3(),
          new TestProblem4(), new TestProblem5(), new TestProblem6()
      }) {

      double previousValueError = Double.NaN;
      double previousTimeError = Double.NaN;
      for (int i = 4; i < 10; ++i) {

        double step = (pb.getFinalTime() - pb.getInitialTime()) * FastMath.pow(2.0, -i);
        FirstOrderIntegrator integ = new MidpointIntegrator(step);
        TestProblemHandler handler = new TestProblemHandler(pb, integ);
        integ.addStepHandler(handler);
        EventHandler[] functions = pb.getEventsHandlers();
        for (int l = 0; l < functions.length; ++l) {
          integ.addEventHandler(functions[l],
                                     Double.POSITIVE_INFINITY, 1.0e-6 * step, 1000);
        }
        double stopTime = integ.integrate(pb,
                                          pb.getInitialTime(), pb.getInitialState(),
                                          pb.getFinalTime(), new double[pb.getDimension()]);
        if (functions.length == 0) {
            Assert.assertEquals(pb.getFinalTime(), stopTime, 1.0e-10);
        }

        double valueError = handler.getMaximalValueError();
        if (i > 4) {
          Assert.assertTrue(valueError < FastMath.abs(previousValueError));
        }
        previousValueError = valueError;

        double timeError = handler.getMaximalTimeError();
        if (i > 4) {
          Assert.assertTrue(timeError <= FastMath.abs(previousTimeError));
        }
        previousTimeError = timeError;

      }

    }

  }

[2338, 'src/test/java', 'org.apache.commons.math3.geometry.euclidean.threed', 'FieldRotationDSTest', 'testRevertVectorOperator', 242, 312]

[2338, 'src/test/java', 'org.apache.commons.math3.geometry.euclidean.threed', 'FieldRotationDSTest', 'testRevertFrameTransform', 314, 384]

    @Test
    public void testRevertVectorOperator() {
        double a = 0.001;
        double b = 0.36;
        double c = 0.48;
        double d = 0.8;
        FieldRotation<DerivativeStructure> r = createRotation(a, b, c, d, true);
        double a2 = a * a;
        double b2 = b * b;
        double c2 = c * c;
        double d2 = d * d;
        double den = (a2 + b2 + c2 + d2) * FastMath.sqrt(a2 + b2 + c2 + d2);
        Assert.assertEquals((b2 + c2 + d2) / den, r.getQ0().getPartialDerivative(1, 0, 0, 0), 1.0e-15);
        Assert.assertEquals(-a * b / den, r.getQ0().getPartialDerivative(0, 1, 0, 0), 1.0e-15);
        Assert.assertEquals(-a * c / den, r.getQ0().getPartialDerivative(0, 0, 1, 0), 1.0e-15);
        Assert.assertEquals(-a * d / den, r.getQ0().getPartialDerivative(0, 0, 0, 1), 1.0e-15);
        Assert.assertEquals(-b * a / den, r.getQ1().getPartialDerivative(1, 0, 0, 0), 1.0e-15);
        Assert.assertEquals((a2 + c2 + d2) / den, r.getQ1().getPartialDerivative(0, 1, 0, 0), 1.0e-15);
        Assert.assertEquals(-b * c / den, r.getQ1().getPartialDerivative(0, 0, 1, 0), 1.0e-15);
        Assert.assertEquals(-b * d / den, r.getQ1().getPartialDerivative(0, 0, 0, 1), 1.0e-15);
        Assert.assertEquals(-c * a / den, r.getQ2().getPartialDerivative(1, 0, 0, 0), 1.0e-15);
        Assert.assertEquals(-c * b / den, r.getQ2().getPartialDerivative(0, 1, 0, 0), 1.0e-15);
        Assert.assertEquals((a2 + b2 + d2) / den, r.getQ2().getPartialDerivative(0, 0, 1, 0), 1.0e-15);
        Assert.assertEquals(-c * d / den, r.getQ2().getPartialDerivative(0, 0, 0, 1), 1.0e-15);
        Assert.assertEquals(-d * a / den, r.getQ3().getPartialDerivative(1, 0, 0, 0), 1.0e-15);
        Assert.assertEquals(-d * b / den, r.getQ3().getPartialDerivative(0, 1, 0, 0), 1.0e-15);
        Assert.assertEquals(-d * c / den, r.getQ3().getPartialDerivative(0, 0, 1, 0), 1.0e-15);
        Assert.assertEquals((a2 + b2 + c2) / den, r.getQ3().getPartialDerivative(0, 0, 0, 1), 1.0e-15);
        FieldRotation<DerivativeStructure> reverted = r.revert();
        FieldRotation<DerivativeStructure> rrT = r.compose(reverted, RotationConvention.VECTOR_OPERATOR);
        checkRotationDS(rrT, 1, 0, 0, 0);
        Assert.assertEquals(0, rrT.getQ0().getPartialDerivative(1, 0, 0, 0), 1.0e-15);
        Assert.assertEquals(0, rrT.getQ0().getPartialDerivative(0, 1, 0, 0), 1.0e-15);
        Assert.assertEquals(0, rrT.getQ0().getPartialDerivative(0, 0, 1, 0), 1.0e-15);
        Assert.assertEquals(0, rrT.getQ0().getPartialDerivative(0, 0, 0, 1), 1.0e-15);
        Assert.assertEquals(0, rrT.getQ1().getPartialDerivative(1, 0, 0, 0), 1.0e-15);
        Assert.assertEquals(0, rrT.getQ1().getPartialDerivative(0, 1, 0, 0), 1.0e-15);
        Assert.assertEquals(0, rrT.getQ1().getPartialDerivative(0, 0, 1, 0), 1.0e-15);
        Assert.assertEquals(0, rrT.getQ1().getPartialDerivative(0, 0, 0, 1), 1.0e-15);
        Assert.assertEquals(0, rrT.getQ2().getPartialDerivative(1, 0, 0, 0), 1.0e-15);
        Assert.assertEquals(0, rrT.getQ2().getPartialDerivative(0, 1, 0, 0), 1.0e-15);
        Assert.assertEquals(0, rrT.getQ2().getPartialDerivative(0, 0, 1, 0), 1.0e-15);
        Assert.assertEquals(0, rrT.getQ2().getPartialDerivative(0, 0, 0, 1), 1.0e-15);
        Assert.assertEquals(0, rrT.getQ3().getPartialDerivative(1, 0, 0, 0), 1.0e-15);
        Assert.assertEquals(0, rrT.getQ3().getPartialDerivative(0, 1, 0, 0), 1.0e-15);
        Assert.assertEquals(0, rrT.getQ3().getPartialDerivative(0, 0, 1, 0), 1.0e-15);
        Assert.assertEquals(0, rrT.getQ3().getPartialDerivative(0, 0, 0, 1), 1.0e-15);
        FieldRotation<DerivativeStructure> rTr = reverted.compose(r, RotationConvention.VECTOR_OPERATOR);
        checkRotationDS(rTr, 1, 0, 0, 0);
        Assert.assertEquals(0, rTr.getQ0().getPartialDerivative(1, 0, 0, 0), 1.0e-15);
        Assert.assertEquals(0, rTr.getQ0().getPartialDerivative(0, 1, 0, 0), 1.0e-15);
        Assert.assertEquals(0, rTr.getQ0().getPartialDerivative(0, 0, 1, 0), 1.0e-15);
        Assert.assertEquals(0, rTr.getQ0().getPartialDerivative(0, 0, 0, 1), 1.0e-15);
        Assert.assertEquals(0, rTr.getQ1().getPartialDerivative(1, 0, 0, 0), 1.0e-15);
        Assert.assertEquals(0, rTr.getQ1().getPartialDerivative(0, 1, 0, 0), 1.0e-15);
        Assert.assertEquals(0, rTr.getQ1().getPartialDerivative(0, 0, 1, 0), 1.0e-15);
        Assert.assertEquals(0, rTr.getQ1().getPartialDerivative(0, 0, 0, 1), 1.0e-15);
        Assert.assertEquals(0, rTr.getQ2().getPartialDerivative(1, 0, 0, 0), 1.0e-15);
        Assert.assertEquals(0, rTr.getQ2().getPartialDerivative(0, 1, 0, 0), 1.0e-15);
        Assert.assertEquals(0, rTr.getQ2().getPartialDerivative(0, 0, 1, 0), 1.0e-15);
        Assert.assertEquals(0, rTr.getQ2().getPartialDerivative(0, 0, 0, 1), 1.0e-15);
        Assert.assertEquals(0, rTr.getQ3().getPartialDerivative(1, 0, 0, 0), 1.0e-15);
        Assert.assertEquals(0, rTr.getQ3().getPartialDerivative(0, 1, 0, 0), 1.0e-15);
        Assert.assertEquals(0, rTr.getQ3().getPartialDerivative(0, 0, 1, 0), 1.0e-15);
        Assert.assertEquals(0, rTr.getQ3().getPartialDerivative(0, 0, 0, 1), 1.0e-15);
        Assert.assertEquals(r.getAngle().getReal(), reverted.getAngle().getReal(), 1.0e-15);
        Assert.assertEquals(-1,
                            FieldVector3D.dotProduct(r.getAxis(RotationConvention.VECTOR_OPERATOR),
                                                     reverted.getAxis(RotationConvention.VECTOR_OPERATOR)).getReal(),
                            1.0e-15);
    }

    @Test
    public void testRevertFrameTransform() {
        double a = 0.001;
        double b = 0.36;
        double c = 0.48;
        double d = 0.8;
        FieldRotation<DerivativeStructure> r = createRotation(a, b, c, d, true);
        double a2 = a * a;
        double b2 = b * b;
        double c2 = c * c;
        double d2 = d * d;
        double den = (a2 + b2 + c2 + d2) * FastMath.sqrt(a2 + b2 + c2 + d2);
        Assert.assertEquals((b2 + c2 + d2) / den, r.getQ0().getPartialDerivative(1, 0, 0, 0), 1.0e-15);
        Assert.assertEquals(-a * b / den, r.getQ0().getPartialDerivative(0, 1, 0, 0), 1.0e-15);
        Assert.assertEquals(-a * c / den, r.getQ0().getPartialDerivative(0, 0, 1, 0), 1.0e-15);
        Assert.assertEquals(-a * d / den, r.getQ0().getPartialDerivative(0, 0, 0, 1), 1.0e-15);
        Assert.assertEquals(-b * a / den, r.getQ1().getPartialDerivative(1, 0, 0, 0), 1.0e-15);
        Assert.assertEquals((a2 + c2 + d2) / den, r.getQ1().getPartialDerivative(0, 1, 0, 0), 1.0e-15);
        Assert.assertEquals(-b * c / den, r.getQ1().getPartialDerivative(0, 0, 1, 0), 1.0e-15);
        Assert.assertEquals(-b * d / den, r.getQ1().getPartialDerivative(0, 0, 0, 1), 1.0e-15);
        Assert.assertEquals(-c * a / den, r.getQ2().getPartialDerivative(1, 0, 0, 0), 1.0e-15);
        Assert.assertEquals(-c * b / den, r.getQ2().getPartialDerivative(0, 1, 0, 0), 1.0e-15);
        Assert.assertEquals((a2 + b2 + d2) / den, r.getQ2().getPartialDerivative(0, 0, 1, 0), 1.0e-15);
        Assert.assertEquals(-c * d / den, r.getQ2().getPartialDerivative(0, 0, 0, 1), 1.0e-15);
        Assert.assertEquals(-d * a / den, r.getQ3().getPartialDerivative(1, 0, 0, 0), 1.0e-15);
        Assert.assertEquals(-d * b / den, r.getQ3().getPartialDerivative(0, 1, 0, 0), 1.0e-15);
        Assert.assertEquals(-d * c / den, r.getQ3().getPartialDerivative(0, 0, 1, 0), 1.0e-15);
        Assert.assertEquals((a2 + b2 + c2) / den, r.getQ3().getPartialDerivative(0, 0, 0, 1), 1.0e-15);
        FieldRotation<DerivativeStructure> reverted = r.revert();
        FieldRotation<DerivativeStructure> rrT = r.compose(reverted, RotationConvention.FRAME_TRANSFORM);
        checkRotationDS(rrT, 1, 0, 0, 0);
        Assert.assertEquals(0, rrT.getQ0().getPartialDerivative(1, 0, 0, 0), 1.0e-15);
        Assert.assertEquals(0, rrT.getQ0().getPartialDerivative(0, 1, 0, 0), 1.0e-15);
        Assert.assertEquals(0, rrT.getQ0().getPartialDerivative(0, 0, 1, 0), 1.0e-15);
        Assert.assertEquals(0, rrT.getQ0().getPartialDerivative(0, 0, 0, 1), 1.0e-15);
        Assert.assertEquals(0, rrT.getQ1().getPartialDerivative(1, 0, 0, 0), 1.0e-15);
        Assert.assertEquals(0, rrT.getQ1().getPartialDerivative(0, 1, 0, 0), 1.0e-15);
        Assert.assertEquals(0, rrT.getQ1().getPartialDerivative(0, 0, 1, 0), 1.0e-15);
        Assert.assertEquals(0, rrT.getQ1().getPartialDerivative(0, 0, 0, 1), 1.0e-15);
        Assert.assertEquals(0, rrT.getQ2().getPartialDerivative(1, 0, 0, 0), 1.0e-15);
        Assert.assertEquals(0, rrT.getQ2().getPartialDerivative(0, 1, 0, 0), 1.0e-15);
        Assert.assertEquals(0, rrT.getQ2().getPartialDerivative(0, 0, 1, 0), 1.0e-15);
        Assert.assertEquals(0, rrT.getQ2().getPartialDerivative(0, 0, 0, 1), 1.0e-15);
        Assert.assertEquals(0, rrT.getQ3().getPartialDerivative(1, 0, 0, 0), 1.0e-15);
        Assert.assertEquals(0, rrT.getQ3().getPartialDerivative(0, 1, 0, 0), 1.0e-15);
        Assert.assertEquals(0, rrT.getQ3().getPartialDerivative(0, 0, 1, 0), 1.0e-15);
        Assert.assertEquals(0, rrT.getQ3().getPartialDerivative(0, 0, 0, 1), 1.0e-15);
        FieldRotation<DerivativeStructure> rTr = reverted.compose(r, RotationConvention.FRAME_TRANSFORM);
        checkRotationDS(rTr, 1, 0, 0, 0);
        Assert.assertEquals(0, rTr.getQ0().getPartialDerivative(1, 0, 0, 0), 1.0e-15);
        Assert.assertEquals(0, rTr.getQ0().getPartialDerivative(0, 1, 0, 0), 1.0e-15);
        Assert.assertEquals(0, rTr.getQ0().getPartialDerivative(0, 0, 1, 0), 1.0e-15);
        Assert.assertEquals(0, rTr.getQ0().getPartialDerivative(0, 0, 0, 1), 1.0e-15);
        Assert.assertEquals(0, rTr.getQ1().getPartialDerivative(1, 0, 0, 0), 1.0e-15);
        Assert.assertEquals(0, rTr.getQ1().getPartialDerivative(0, 1, 0, 0), 1.0e-15);
        Assert.assertEquals(0, rTr.getQ1().getPartialDerivative(0, 0, 1, 0), 1.0e-15);
        Assert.assertEquals(0, rTr.getQ1().getPartialDerivative(0, 0, 0, 1), 1.0e-15);
        Assert.assertEquals(0, rTr.getQ2().getPartialDerivative(1, 0, 0, 0), 1.0e-15);
        Assert.assertEquals(0, rTr.getQ2().getPartialDerivative(0, 1, 0, 0), 1.0e-15);
        Assert.assertEquals(0, rTr.getQ2().getPartialDerivative(0, 0, 1, 0), 1.0e-15);
        Assert.assertEquals(0, rTr.getQ2().getPartialDerivative(0, 0, 0, 1), 1.0e-15);
        Assert.assertEquals(0, rTr.getQ3().getPartialDerivative(1, 0, 0, 0), 1.0e-15);
        Assert.assertEquals(0, rTr.getQ3().getPartialDerivative(0, 1, 0, 0), 1.0e-15);
        Assert.assertEquals(0, rTr.getQ3().getPartialDerivative(0, 0, 1, 0), 1.0e-15);
        Assert.assertEquals(0, rTr.getQ3().getPartialDerivative(0, 0, 0, 1), 1.0e-15);
        Assert.assertEquals(r.getAngle().getReal(), reverted.getAngle().getReal(), 1.0e-15);
        Assert.assertEquals(-1,
                            FieldVector3D.dotProduct(r.getAxis(RotationConvention.FRAME_TRANSFORM),
                                                     reverted.getAxis(RotationConvention.FRAME_TRANSFORM)).getReal(),
                            1.0e-15);
    }

[2345, 'src/test/java', 'org.apache.commons.math3.util', 'MathArraysTest', 'testL1DistanceInt', 140, 145]

[2345, 'src/test/java', 'org.apache.commons.math3.util', 'MathArraysTest', 'testLInfDistanceInt', 168, 173]

    @Test
    public void testL1DistanceInt() {
        int[] p1 = { 3, 0 };
        int[] p2 = { 0, 4 };
        Assert.assertEquals(7, MathArrays.distance1(p1, p2));
    }

    @Test
    public void testLInfDistanceInt() {
        int[] p1 = { 3, 0 };
        int[] p2 = { 0, 4 };
        Assert.assertEquals(4, MathArrays.distanceInf(p1, p2));
    }

[2369, 'src/test/java', 'org.apache.commons.math3.stat', 'CertifiedDataTest', 'testSummaryStatistics', 42, 64]

[2369, 'src/test/java', 'org.apache.commons.math3.stat', 'CertifiedDataTest', 'testDescriptiveStatistics', 70, 94]

    @Test
    public void testSummaryStatistics() throws Exception {
        SummaryStatistics u = new SummaryStatistics();
        loadStats("data/PiDigits.txt", u);
        Assert.assertEquals("PiDigits: std", std, u.getStandardDeviation(), 1E-13);
        Assert.assertEquals("PiDigits: mean", mean, u.getMean(), 1E-13);

        loadStats("data/Mavro.txt", u);
        Assert.assertEquals("Mavro: std", std, u.getStandardDeviation(), 1E-14);
        Assert.assertEquals("Mavro: mean", mean, u.getMean(), 1E-14);

        loadStats("data/Michelso.txt", u);
        Assert.assertEquals("Michelso: std", std, u.getStandardDeviation(), 1E-13);
        Assert.assertEquals("Michelso: mean", mean, u.getMean(), 1E-13);

        loadStats("data/NumAcc1.txt", u);
        Assert.assertEquals("NumAcc1: std", std, u.getStandardDeviation(), 1E-14);
        Assert.assertEquals("NumAcc1: mean", mean, u.getMean(), 1E-14);

        loadStats("data/NumAcc2.txt", u);
        Assert.assertEquals("NumAcc2: std", std, u.getStandardDeviation(), 1E-14);
        Assert.assertEquals("NumAcc2: mean", mean, u.getMean(), 1E-14);
    }

    @Test
    public void testDescriptiveStatistics() throws Exception {

        DescriptiveStatistics u = new DescriptiveStatistics();

        loadStats("data/PiDigits.txt", u);
        Assert.assertEquals("PiDigits: std", std, u.getStandardDeviation(), 1E-14);
        Assert.assertEquals("PiDigits: mean", mean, u.getMean(), 1E-14);

        loadStats("data/Mavro.txt", u);
        Assert.assertEquals("Mavro: std", std, u.getStandardDeviation(), 1E-14);
        Assert.assertEquals("Mavro: mean", mean, u.getMean(), 1E-14);

        loadStats("data/Michelso.txt", u);
        Assert.assertEquals("Michelso: std", std, u.getStandardDeviation(), 1E-14);
        Assert.assertEquals("Michelso: mean", mean, u.getMean(), 1E-14);

        loadStats("data/NumAcc1.txt", u);
        Assert.assertEquals("NumAcc1: std", std, u.getStandardDeviation(), 1E-14);
        Assert.assertEquals("NumAcc1: mean", mean, u.getMean(), 1E-14);

        loadStats("data/NumAcc2.txt", u);
        Assert.assertEquals("NumAcc2: std", std, u.getStandardDeviation(), 1E-14);
        Assert.assertEquals("NumAcc2: mean", mean, u.getMean(), 1E-14);
    }

[2374, 'src/test/java', 'org.apache.commons.math3.fitting', 'GaussianCurveFitterTest', 'testFit02', 235, 238]

[2374, 'src/test/java', 'org.apache.commons.math3.fitting', 'HarmonicCurveFitterTest', 'testPreconditions1', 35, 38]

    @Test(expected=MathIllegalArgumentException.class)
    public void testFit02() {
        GaussianCurveFitter.create().fit(new WeightedObservedPoints().toList());
    }

    @Test(expected=NumberIsTooSmallException.class)
    public void testPreconditions1() {
        HarmonicCurveFitter.create().fit(new WeightedObservedPoints().toList());
    }

[2375, 'src/test/java', 'org.apache.commons.math3.complex', 'ComplexTest', 'testConjugateNaN', 167, 171]

[2375, 'src/test/java', 'org.apache.commons.math3.complex', 'ComplexTest', 'testNegateNaN', 428, 432]

    @Test
    public void testConjugateNaN() {
        Complex z = Complex.NaN.conjugate();
        Assert.assertTrue(z.isNaN());
    }

    @Test
    public void testNegateNaN() {
        Complex z = Complex.NaN.negate();
        Assert.assertTrue(z.isNaN());
    }

[2383, 'src/test/java', 'org.apache.commons.math3.ode.nonstiff', 'ClassicalRungeKuttaStepInterpolatorTest', 'serialization', 51, 96]

[2383, 'src/test/java', 'org.apache.commons.math3.ode.nonstiff', 'ThreeEighthesStepInterpolatorTest', 'serialization', 51, 96]

  @Test
  public void serialization()
    throws IOException, ClassNotFoundException,
           DimensionMismatchException, NumberIsTooSmallException,
           MaxCountExceededException, NoBracketingException  {

    TestProblem3 pb = new TestProblem3(0.9);
    double step = (pb.getFinalTime() - pb.getInitialTime()) * 0.0003;
    ClassicalRungeKuttaIntegrator integ = new ClassicalRungeKuttaIntegrator(step);
    integ.addStepHandler(new ContinuousOutputModel());
    integ.integrate(pb,
                    pb.getInitialTime(), pb.getInitialState(),
                    pb.getFinalTime(), new double[pb.getDimension()]);

    ByteArrayOutputStream bos = new ByteArrayOutputStream();
    ObjectOutputStream    oos = new ObjectOutputStream(bos);
    for (StepHandler handler : integ.getStepHandlers()) {
        oos.writeObject(handler);
    }

    Assert.assertTrue(bos.size () > 880000);
    Assert.assertTrue(bos.size () < 900000);

    ByteArrayInputStream  bis = new ByteArrayInputStream(bos.toByteArray());
    ObjectInputStream     ois = new ObjectInputStream(bis);
    ContinuousOutputModel cm  = (ContinuousOutputModel) ois.readObject();

    Random random = new Random(347588535632l);
    double maxError = 0.0;
    for (int i = 0; i < 1000; ++i) {
      double r = random.nextDouble();
      double time = r * pb.getInitialTime() + (1.0 - r) * pb.getFinalTime();
      cm.setInterpolatedTime(time);
      double[] interpolatedY = cm.getInterpolatedState ();
      double[] theoreticalY  = pb.computeTheoreticalState(time);
      double dx = interpolatedY[0] - theoreticalY[0];
      double dy = interpolatedY[1] - theoreticalY[1];
      double error = dx * dx + dy * dy;
      if (error > maxError) {
        maxError = error;
      }
    }

    Assert.assertTrue(maxError > 0.005);

  }

  @Test
  public void serialization()
    throws IOException, ClassNotFoundException,
           DimensionMismatchException, NumberIsTooSmallException,
           MaxCountExceededException, NoBracketingException {

    TestProblem3 pb = new TestProblem3(0.9);
    double step = (pb.getFinalTime() - pb.getInitialTime()) * 0.0003;
    ThreeEighthesIntegrator integ = new ThreeEighthesIntegrator(step);
    integ.addStepHandler(new ContinuousOutputModel());
    integ.integrate(pb,
                    pb.getInitialTime(), pb.getInitialState(),
                    pb.getFinalTime(), new double[pb.getDimension()]);

    ByteArrayOutputStream bos = new ByteArrayOutputStream();
    ObjectOutputStream    oos = new ObjectOutputStream(bos);
    for (StepHandler handler : integ.getStepHandlers()) {
        oos.writeObject(handler);
    }

    Assert.assertTrue(bos.size () > 880000);
    Assert.assertTrue(bos.size () < 900000);

    ByteArrayInputStream  bis = new ByteArrayInputStream(bos.toByteArray());
    ObjectInputStream     ois = new ObjectInputStream(bis);
    ContinuousOutputModel cm  = (ContinuousOutputModel) ois.readObject();

    Random random = new Random(347588535632l);
    double maxError = 0.0;
    for (int i = 0; i < 1000; ++i) {
      double r = random.nextDouble();
      double time = r * pb.getInitialTime() + (1.0 - r) * pb.getFinalTime();
      cm.setInterpolatedTime(time);
      double[] interpolatedY = cm.getInterpolatedState ();
      double[] theoreticalY  = pb.computeTheoreticalState(time);
      double dx = interpolatedY[0] - theoreticalY[0];
      double dy = interpolatedY[1] - theoreticalY[1];
      double error = dx * dx + dy * dy;
      if (error > maxError) {
        maxError = error;
      }
    }

    Assert.assertTrue(maxError > 0.005);

  }

[2386, 'src/test/java', 'org.apache.commons.math3.linear', 'Array2DRowRealMatrixTest', 'testExamples', 418, 449]

[2386, 'src/test/java', 'org.apache.commons.math3.linear', 'BlockRealMatrixTest', 'testExamples', 463, 494]

    @Test
    public void testExamples() {
        // Create a real matrix with two rows and three columns
        double[][] matrixData = { {1d,2d,3d}, {2d,5d,3d}};
        RealMatrix m = new Array2DRowRealMatrix(matrixData);
        // One more with three rows, two columns
        double[][] matrixData2 = { {1d,2d}, {2d,5d}, {1d, 7d}};
        RealMatrix n = new Array2DRowRealMatrix(matrixData2);
        // Now multiply m by n
        RealMatrix p = m.multiply(n);
        Assert.assertEquals(2, p.getRowDimension());
        Assert.assertEquals(2, p.getColumnDimension());
        // Invert p
        RealMatrix pInverse = new LUDecomposition(p).getSolver().getInverse();
        Assert.assertEquals(2, pInverse.getRowDimension());
        Assert.assertEquals(2, pInverse.getColumnDimension());

        // Solve example
        double[][] coefficientsData = {{2, 3, -2}, {-1, 7, 6}, {4, -3, -5}};
        RealMatrix coefficients = new Array2DRowRealMatrix(coefficientsData);
        RealVector constants = new ArrayRealVector(new double[]{1, -2, 1}, false);
        RealVector solution = new LUDecomposition(coefficients).getSolver().solve(constants);
        final double cst0 = constants.getEntry(0);
        final double cst1 = constants.getEntry(1);
        final double cst2 = constants.getEntry(2);
        final double sol0 = solution.getEntry(0);
        final double sol1 = solution.getEntry(1);
        final double sol2 = solution.getEntry(2);
        Assert.assertEquals(2 * sol0 + 3 * sol1 -2 * sol2, cst0, 1E-12);
        Assert.assertEquals(-1 * sol0 + 7 * sol1 + 6 * sol2, cst1, 1E-12);
        Assert.assertEquals(4 * sol0 - 3 * sol1 -5 * sol2, cst2, 1E-12);
    }

    @Test
    public void testExamples() {
        // Create a real matrix with two rows and three columns
        double[][] matrixData = { {1d,2d,3d}, {2d,5d,3d}};
        RealMatrix m = new BlockRealMatrix(matrixData);
        // One more with three rows, two columns
        double[][] matrixData2 = { {1d,2d}, {2d,5d}, {1d, 7d}};
        RealMatrix n = new BlockRealMatrix(matrixData2);
        // Now multiply m by n
        RealMatrix p = m.multiply(n);
        Assert.assertEquals(2, p.getRowDimension());
        Assert.assertEquals(2, p.getColumnDimension());
        // Invert p
        RealMatrix pInverse = new LUDecomposition(p).getSolver().getInverse();
        Assert.assertEquals(2, pInverse.getRowDimension());
        Assert.assertEquals(2, pInverse.getColumnDimension());

        // Solve example
        double[][] coefficientsData = {{2, 3, -2}, {-1, 7, 6}, {4, -3, -5}};
        RealMatrix coefficients = new BlockRealMatrix(coefficientsData);
        RealVector constants = new ArrayRealVector(new double[]{1, -2, 1}, false);
        RealVector solution = new LUDecomposition(coefficients).getSolver().solve(constants);
        final double cst0 = constants.getEntry(0);
        final double cst1 = constants.getEntry(1);
        final double cst2 = constants.getEntry(2);
        final double sol0 = solution.getEntry(0);
        final double sol1 = solution.getEntry(1);
        final double sol2 = solution.getEntry(2);
        Assert.assertEquals(2 * sol0 + 3 * sol1 -2 * sol2, cst0, 1E-12);
        Assert.assertEquals(-1 * sol0 + 7 * sol1 + 6 * sol2, cst1, 1E-12);
        Assert.assertEquals(4 * sol0 - 3 * sol1 -5 * sol2, cst2, 1E-12);
    }

[2393, 'src/test/java', 'org.apache.commons.math3.stat', 'FrequencyTest', 'testNonComparableCumPct', 201, 206]

[2393, 'src/test/java', 'org.apache.commons.math3.stat', 'FrequencyTest', 'testNonComparablePct', 208, 213]

    @Test
    public void testNonComparableCumPct() {
        f.addValue("a");
        Assert.assertEquals("cum freq, single entry", 1.0d, f.getCumPct("a"),TOLERANCE);
        Assert.assertEquals("cum freq, single entry non comparable", 0.0d, f.getCumPct(100),TOLERANCE);
    }

    @Test
    public void testNonComparablePct() {
        f.addValue("a");
        Assert.assertEquals("cum freq, single entry", 1.0d, f.getPct("a"),TOLERANCE);
        Assert.assertEquals("cum freq, single entry non comparable", 0.0d, f.getPct(100),TOLERANCE);
    }

[2395, 'src/test/java', 'org.apache.commons.math3.complex', 'ComplexTest', 'testNthRoot_cornercase_thirdRoot_imaginaryPartEmpty', 1279, 1297]

[2395, 'src/test/java', 'org.apache.commons.math3.complex', 'ComplexTest', 'testNthRoot_cornercase_thirdRoot_realPartZero', 1311, 1328]

    @Test
    public void testNthRoot_cornercase_thirdRoot_imaginaryPartEmpty() {
        // The number 8 has three third roots. One we all already know is the number 2.
        // But there are two more complex roots.
        Complex z = new Complex(8,0);
        // The List holding all third roots
        Complex[] thirdRootsOfZ = z.nthRoot(3).toArray(new Complex[0]);
        // Returned Collection must not be empty!
        Assert.assertEquals(3, thirdRootsOfZ.length);
        // test z_0
        Assert.assertEquals(2.0,                thirdRootsOfZ[0].getReal(),      1.0e-5);
        Assert.assertEquals(0.0,                thirdRootsOfZ[0].getImaginary(), 1.0e-5);
        // test z_1
        Assert.assertEquals(-1.0,               thirdRootsOfZ[1].getReal(),      1.0e-5);
        Assert.assertEquals(1.7320508075688774, thirdRootsOfZ[1].getImaginary(), 1.0e-5);
        // test z_2
        Assert.assertEquals(-1.0,               thirdRootsOfZ[2].getReal(),      1.0e-5);
        Assert.assertEquals(-1.732050807568877, thirdRootsOfZ[2].getImaginary(), 1.0e-5);
    }

    @Test
    public void testNthRoot_cornercase_thirdRoot_realPartZero() {
        // complex number with only imaginary part
        Complex z = new Complex(0,2);
        // The List holding all third roots
        Complex[] thirdRootsOfZ = z.nthRoot(3).toArray(new Complex[0]);
        // Returned Collection must not be empty!
        Assert.assertEquals(3, thirdRootsOfZ.length);
        // test z_0
        Assert.assertEquals(1.0911236359717216,      thirdRootsOfZ[0].getReal(),      1.0e-5);
        Assert.assertEquals(0.6299605249474365,      thirdRootsOfZ[0].getImaginary(), 1.0e-5);
        // test z_1
        Assert.assertEquals(-1.0911236359717216,     thirdRootsOfZ[1].getReal(),      1.0e-5);
        Assert.assertEquals(0.6299605249474365,      thirdRootsOfZ[1].getImaginary(), 1.0e-5);
        // test z_2
        Assert.assertEquals(-2.3144374213981936E-16, thirdRootsOfZ[2].getReal(),      1.0e-5);
        Assert.assertEquals(-1.2599210498948732,     thirdRootsOfZ[2].getImaginary(), 1.0e-5);
    }

[2399, 'src/test/java', 'org.apache.commons.math3.analysis.interpolation', 'BicubicInterpolatingFunctionTest', 'testIsValidPoint', 115, 173]

[2399, 'src/test/java', 'org.apache.commons.math3.analysis.interpolation', 'BicubicSplineInterpolatingFunctionTest', 'testIsValidPoint', 612, 670]

    @Test
    public void testIsValidPoint() {
        final double xMin = -12;
        final double xMax = 34;
        final double yMin = 5;
        final double yMax = 67;
        final double[] xval = new double[] { xMin, xMax };
        final double[] yval = new double[] { yMin, yMax };
        final double[][] f = new double[][] { { 1, 2 },
                                              { 3, 4 } };
        final double[][] dFdX = f;
        final double[][] dFdY = f;
        final double[][] dFdXdY = f;

        final BicubicInterpolatingFunction bcf
            = new BicubicInterpolatingFunction(xval, yval, f,
                                                     dFdX, dFdY, dFdXdY);

        double x, y;

        x = xMin;
        y = yMin;
        Assert.assertTrue(bcf.isValidPoint(x, y));
        // Ensure that no exception is thrown.
        bcf.value(x, y);

        x = xMax;
        y = yMax;
        Assert.assertTrue(bcf.isValidPoint(x, y));
        // Ensure that no exception is thrown.
        bcf.value(x, y);

        final double xRange = xMax - xMin;
        final double yRange = yMax - yMin;
        x = xMin + xRange / 3.4;
        y = yMin + yRange / 1.2;
        Assert.assertTrue(bcf.isValidPoint(x, y));
        // Ensure that no exception is thrown.
        bcf.value(x, y);

        final double small = 1e-8;
        x = xMin - small;
        y = yMax;
        Assert.assertFalse(bcf.isValidPoint(x, y));
        // Ensure that an exception would have been thrown.
        try {
            bcf.value(x, y);
            Assert.fail("OutOfRangeException expected");
        } catch (OutOfRangeException expected) {}

        x = xMin;
        y = yMax + small;
        Assert.assertFalse(bcf.isValidPoint(x, y));
        // Ensure that an exception would have been thrown.
        try {
            bcf.value(x, y);
            Assert.fail("OutOfRangeException expected");
        } catch (OutOfRangeException expected) {}
    }

    @Test
    public void testIsValidPoint() {
        final double xMin = -12;
        final double xMax = 34;
        final double yMin = 5;
        final double yMax = 67;
        final double[] xval = new double[] { xMin, xMax };
        final double[] yval = new double[] { yMin, yMax };
        final double[][] f = new double[][] { { 1, 2 },
                                              { 3, 4 } };
        final double[][] dFdX = f;
        final double[][] dFdY = f;
        final double[][] dFdXdY = f;

        final BicubicSplineInterpolatingFunction bcf
            = new BicubicSplineInterpolatingFunction(xval, yval, f,
                                                     dFdX, dFdY, dFdXdY);

        double x, y;

        x = xMin;
        y = yMin;
        Assert.assertTrue(bcf.isValidPoint(x, y));
        // Ensure that no exception is thrown.
        bcf.value(x, y);

        x = xMax;
        y = yMax;
        Assert.assertTrue(bcf.isValidPoint(x, y));
        // Ensure that no exception is thrown.
        bcf.value(x, y);

        final double xRange = xMax - xMin;
        final double yRange = yMax - yMin;
        x = xMin + xRange / 3.4;
        y = yMin + yRange / 1.2;
        Assert.assertTrue(bcf.isValidPoint(x, y));
        // Ensure that no exception is thrown.
        bcf.value(x, y);

        final double small = 1e-8;
        x = xMin - small;
        y = yMax;
        Assert.assertFalse(bcf.isValidPoint(x, y));
        // Ensure that an exception would have been thrown.
        try {
            bcf.value(x, y);
            Assert.fail("OutOfRangeException expected");
        } catch (OutOfRangeException expected) {}

        x = xMin;
        y = yMax + small;
        Assert.assertFalse(bcf.isValidPoint(x, y));
        // Ensure that an exception would have been thrown.
        try {
            bcf.value(x, y);
            Assert.fail("OutOfRangeException expected");
        } catch (OutOfRangeException expected) {}
    }

[2404, 'src/test/java', 'org.apache.commons.math3.linear', 'BlockFieldMatrixTest', 'testSetColumnVector', 952, 971]

[2404, 'src/test/java', 'org.apache.commons.math3.linear', 'FieldMatrixImplTest', 'testSetColumnVector', 732, 751]

    @Test
    public void testSetColumnVector() {
        FieldMatrix<Fraction> m = new BlockFieldMatrix<Fraction>(subTestData);
        FieldVector<Fraction> mColumn3 = columnToVector(subColumn3);
        Assert.assertNotSame(mColumn3, m.getColumnVector(1));
        m.setColumnVector(1, mColumn3);
        Assert.assertEquals(mColumn3, m.getColumnVector(1));
        try {
            m.setColumnVector(-1, mColumn3);
            Assert.fail("Expecting OutOfRangeException");
        } catch (OutOfRangeException ex) {
            // expected
        }
        try {
            m.setColumnVector(0, new ArrayFieldVector<Fraction>(FractionField.getInstance(), 5));
            Assert.fail("Expecting MatrixDimensionMismatchException");
        } catch (MatrixDimensionMismatchException ex) {
            // expected
        }
    }

    @Test
    public void testSetColumnVector() {
        FieldMatrix<Fraction> m = new Array2DRowFieldMatrix<Fraction>(subTestData);
        FieldVector<Fraction> mColumn3 = columnToVector(subColumn3);
        Assert.assertNotSame(mColumn3, m.getColumnVector(1));
        m.setColumnVector(1, mColumn3);
        Assert.assertEquals(mColumn3, m.getColumnVector(1));
        try {
            m.setColumnVector(-1, mColumn3);
            Assert.fail("Expecting OutOfRangeException");
        } catch (OutOfRangeException ex) {
            // expected
        }
        try {
            m.setColumnVector(0, new ArrayFieldVector<Fraction>(FractionField.getInstance(), 5));
            Assert.fail("Expecting MatrixDimensionMismatchException");
        } catch (MatrixDimensionMismatchException ex) {
            // expected
        }
    }

[2443, 'src/test/java', 'org.apache.commons.math3.geometry.euclidean.threed', 'FieldRotationDfpTest', 'testRevertVectorOperator', 188, 205]

[2443, 'src/test/java', 'org.apache.commons.math3.geometry.euclidean.threed', 'FieldRotationDfpTest', 'testRevertFrameTransform', 207, 224]

    @Test
    public void testRevertVectorOperator() {
        double a = 0.001;
        double b = 0.36;
        double c = 0.48;
        double d = 0.8;
        FieldRotation<Dfp> r = createRotation(a, b, c, d, true);
        FieldRotation<Dfp> reverted = r.revert();
        FieldRotation<Dfp> rrT = r.compose(reverted, RotationConvention.VECTOR_OPERATOR);
        checkRotationDS(rrT, 1, 0, 0, 0);
        FieldRotation<Dfp> rTr = reverted.compose(r, RotationConvention.VECTOR_OPERATOR);
        checkRotationDS(rTr, 1, 0, 0, 0);
        Assert.assertEquals(r.getAngle().getReal(), reverted.getAngle().getReal(), 1.0e-15);
        Assert.assertEquals(-1,
                            FieldVector3D.dotProduct(r.getAxis(RotationConvention.VECTOR_OPERATOR),
                                                     reverted.getAxis(RotationConvention.VECTOR_OPERATOR)).getReal(),
                            1.0e-15);
    }

    @Test
    public void testRevertFrameTransform() {
        double a = 0.001;
        double b = 0.36;
        double c = 0.48;
        double d = 0.8;
        FieldRotation<Dfp> r = createRotation(a, b, c, d, true);
        FieldRotation<Dfp> reverted = r.revert();
        FieldRotation<Dfp> rrT = r.compose(reverted, RotationConvention.FRAME_TRANSFORM);
        checkRotationDS(rrT, 1, 0, 0, 0);
        FieldRotation<Dfp> rTr = reverted.compose(r, RotationConvention.FRAME_TRANSFORM);
        checkRotationDS(rTr, 1, 0, 0, 0);
        Assert.assertEquals(r.getAngle().getReal(), reverted.getAngle().getReal(), 1.0e-15);
        Assert.assertEquals(-1,
                            FieldVector3D.dotProduct(r.getAxis(RotationConvention.FRAME_TRANSFORM),
                                                     reverted.getAxis(RotationConvention.FRAME_TRANSFORM)).getReal(),
                            1.0e-15);
    }

[2460, 'src/test/java', 'org.apache.commons.math3.stat.correlation', 'KendallsCorrelationTest', 'testAllTiesInX', 181, 190]

[2460, 'src/test/java', 'org.apache.commons.math3.stat.correlation', 'KendallsCorrelationTest', 'testAllTiesInY', 192, 201]

    @Test
    public void testAllTiesInX() {
        final int length = 10;
        final double[] xArray = new double[length];
        final double[] yArray = new double[length];
        for (int i = 0; i < length; i++) {
            xArray[i] = i;
        }
        Assert.assertEquals(Double.NaN, correlation.correlation(xArray, yArray), 0);
    }

    @Test
    public void testAllTiesInY() {
        final int length = 10;
        final double[] xArray = new double[length];
        final double[] yArray = new double[length];
        for (int i = 0; i < length; i++) {
            yArray[i] = i;
        }
        Assert.assertEquals(Double.NaN, correlation.correlation(xArray, yArray), 0);
    }

[2469, 'src/test/java', 'org.apache.commons.math3.distribution', 'LaplaceDistributionTest', 'testParameters', 28, 33]

[2469, 'src/test/java', 'org.apache.commons.math3.distribution', 'LogisticsDistributionTest', 'testParameters', 28, 33]

    @Test
    public void testParameters() {
        LaplaceDistribution d = makeDistribution();
        Assert.assertEquals(0, d.getLocation(), Precision.EPSILON);
        Assert.assertEquals(1, d.getScale(), Precision.EPSILON);
    }

    @Test
    public void testParameters() {
        LogisticDistribution d = makeDistribution();
        Assert.assertEquals(2, d.getLocation(), Precision.EPSILON);
        Assert.assertEquals(5, d.getScale(), Precision.EPSILON);
    }

[2470, 'src/test/java', 'org.apache.commons.math3.distribution', 'GumbelDistributionTest', 'testParameters', 28, 33]

[2470, 'src/test/java', 'org.apache.commons.math3.distribution', 'NakagamiDistributionTest', 'testParameters', 28, 33]

    @Test
    public void testParameters() {
        GumbelDistribution d = makeDistribution();
        Assert.assertEquals(0.5, d.getLocation(), Precision.EPSILON);
        Assert.assertEquals(2, d.getScale(), Precision.EPSILON);
    }

    @Test
    public void testParameters() {
        NakagamiDistribution d = makeDistribution();
        Assert.assertEquals(0.5, d.getShape(), Precision.EPSILON);
        Assert.assertEquals(1, d.getScale(), Precision.EPSILON);
    }

[2474, 'src/test/java', 'org.apache.commons.math3.linear', 'BlockFieldMatrixTest', 'testScalarAdd', 367, 372]

[2474, 'src/test/java', 'org.apache.commons.math3.linear', 'FieldMatrixImplTest', 'testScalarAdd', 261, 265]

    @Test
    public void testScalarAdd() {
        FieldMatrix<Fraction> m = new BlockFieldMatrix<Fraction>(testData);
        TestUtils.assertEquals(new BlockFieldMatrix<Fraction>(testDataPlus2),
                               m.scalarAdd(new Fraction(2)));
    }

    @Test
    public void testScalarAdd() {
        FieldMatrix<Fraction> m = new Array2DRowFieldMatrix<Fraction>(testData);
        TestUtils.assertEquals(new Array2DRowFieldMatrix<Fraction>(testDataPlus2), m.scalarAdd(new Fraction(2)));
    }

[2491, 'src/test/java', 'org.apache.commons.math3.linear', 'SparseFieldMatrixTest', 'testOperate', 261, 275]

[2491, 'src/test/java', 'org.apache.commons.math3.linear', 'SparseFieldMatrixTest', 'testPremultiplyVector', 302, 316]

    @Test
    public void testOperate() {
        FieldMatrix<Fraction> m = createSparseMatrix(id);
        assertClose("identity operate", testVector, m.operate(testVector),
                entryTolerance);
        assertClose("identity operate", testVector, m.operate(
                new ArrayFieldVector<Fraction>(testVector)).toArray(), entryTolerance);
        m = createSparseMatrix(bigSingular);
        try {
            m.operate(testVector);
            Assert.fail("Expecting illegalArgumentException");
        } catch (MathIllegalArgumentException ex) {
            // ignored
        }
    }

    @Test
    public void testPremultiplyVector() {
        FieldMatrix<Fraction> m = createSparseMatrix(testData);
        assertClose("premultiply", m.preMultiply(testVector), preMultTest,
            normTolerance);
        assertClose("premultiply", m.preMultiply(
            new ArrayFieldVector<Fraction>(testVector).getData()), preMultTest, normTolerance);
        m = createSparseMatrix(bigSingular);
        try {
            m.preMultiply(testVector);
            Assert.fail("expecting MathIllegalArgumentException");
        } catch (MathIllegalArgumentException ex) {
            // ignored
        }
    }

[2520, 'src/test/java', 'org.apache.commons.math3.ode.events', 'EventFilterTest', 'testIncreasingOnly', 132, 157]

[2520, 'src/test/java', 'org.apache.commons.math3.ode.events', 'EventFilterTest', 'testDecreasingOnly', 159, 184]

    @Test
    public void testIncreasingOnly()
        throws DimensionMismatchException, NumberIsTooSmallException,
               MaxCountExceededException, NoBracketingException {
        double e = 1e-15;
        FirstOrderIntegrator integrator;
        integrator = new DormandPrince853Integrator(1.0e-3, 100.0, 1e-7, 1e-7);
        Event allEvents = new Event(true, true);
        integrator.addEventHandler(allEvents, 0.1, e, 1000,
                                   new BracketingNthOrderBrentSolver(1.0e-7, 5));
        Event onlyIncreasing = new Event(false, true);
        integrator.addEventHandler(new EventFilter(onlyIncreasing,
                                                   FilterType.TRIGGER_ONLY_INCREASING_EVENTS),
                                   0.1, e, 100,
                                   new BracketingNthOrderBrentSolver(1.0e-7, 5));
        double t0 = 0.5 * FastMath.PI;
        double tEnd = 5.5 * FastMath.PI;
        double[] y = { 0.0, 1.0 };
        Assert.assertEquals(tEnd,
                            integrator.integrate(new SineCosine(), t0, y, tEnd, y),
                            1.0e-7);

        Assert.assertEquals(5, allEvents.getEventCount());
        Assert.assertEquals(2, onlyIncreasing.getEventCount());

    }

    @Test
    public void testDecreasingOnly()
        throws DimensionMismatchException, NumberIsTooSmallException,
               MaxCountExceededException, NoBracketingException {
        double e = 1e-15;
        FirstOrderIntegrator integrator;
        integrator = new DormandPrince853Integrator(1.0e-3, 100.0, 1e-7, 1e-7);
        Event allEvents = new Event(true, true);
        integrator.addEventHandler(allEvents, 0.1, e, 1000,
                                   new BracketingNthOrderBrentSolver(1.0e-7, 5));
        Event onlyDecreasing = new Event(true, false);
        integrator.addEventHandler(new EventFilter(onlyDecreasing,
                                                   FilterType.TRIGGER_ONLY_DECREASING_EVENTS),
                                   0.1, e, 1000,
                                   new BracketingNthOrderBrentSolver(1.0e-7, 5));
        double t0 = 0.5 * FastMath.PI;
        double tEnd = 5.5 * FastMath.PI;
        double[] y = { 0.0, 1.0 };
        Assert.assertEquals(tEnd,
                            integrator.integrate(new SineCosine(), t0, y, tEnd, y),
                            1.0e-7);

        Assert.assertEquals(5, allEvents.getEventCount());
        Assert.assertEquals(3, onlyDecreasing.getEventCount());

    }

[2531, 'src/test/java', 'org.apache.commons.math3.linear', 'SingularValueDecompositionTest', 'testStability1', 252, 261]

[2531, 'src/test/java', 'org.apache.commons.math3.linear', 'SingularValueDecompositionTest', 'testStability2', 264, 273]

    @Test
    public void testStability1() {
        RealMatrix m = new Array2DRowRealMatrix(201, 201);
        loadRealMatrix(m,"matrix1.csv");
        try {
            new SingularValueDecomposition(m);
        } catch (Exception e) {
            Assert.fail("Exception whilst constructing SVD");
        }
    }

    @Test
    public void testStability2() {
        RealMatrix m = new Array2DRowRealMatrix(7, 168);
        loadRealMatrix(m,"matrix2.csv");
        try {
            new SingularValueDecomposition(m);
        } catch (Throwable e) {
            Assert.fail("Exception whilst constructing SVD");
        }
    }

[2534, 'src/test/java', 'org.apache.commons.math3.geometry.euclidean.oned', 'Vector1DFormatAbstractTest', 'testSimpleNoDecimals', 44, 50]

[2534, 'src/test/java', 'org.apache.commons.math3.geometry.euclidean.oned', 'Vector1DFormatAbstractTest', 'testNonDefaultSetting', 82, 88]

    @Test
    public void testSimpleNoDecimals() {
        Vector1D c = new Vector1D(1);
        String expected = "{1}";
        String actual = vector1DFormat.format(c);
        Assert.assertEquals(expected, actual);
    }

    @Test
    public void testNonDefaultSetting() {
        Vector1D c = new Vector1D(1);
        String expected = "[1]";
        String actual = vector1DFormatSquare.format(c);
        Assert.assertEquals(expected, actual);
    }

[2540, 'src/test/java', 'org.apache.commons.math3.stat.correlation', 'StorelessCovarianceTest', 'testLonglySimpleVar', 116, 124]

[2540, 'src/test/java', 'org.apache.commons.math3.stat.correlation', 'StorelessCovarianceTest', 'testLonglySimpleCov', 126, 134]

    @Test
    public void testLonglySimpleVar(){
        double rCov = 12333921.73333333246;
        StorelessBivariateCovariance cov = new StorelessBivariateCovariance();
        for(int i=0;i<longleyDataSimple.length;i++){
            cov.increment(longleyDataSimple[i][0],longleyDataSimple[i][0]);
        }
        TestUtils.assertEquals("simple covariance test", rCov, cov.getResult(), 10E-7);
    }

    @Test
    public void testLonglySimpleCov(){
        double rCov = 36796.660000;
        StorelessBivariateCovariance cov = new StorelessBivariateCovariance();
        for(int i=0;i<longleyDataSimple.length;i++){
            cov.increment(longleyDataSimple[i][0], longleyDataSimple[i][1]);
        }
        TestUtils.assertEquals("simple covariance test", rCov, cov.getResult(), 10E-7);
    }

[2545, 'src/test/java', 'org.apache.commons.math3.linear', 'RealMatrixFormatAbstractTest', 'testParseManyComponents', 326, 330]

[2545, 'src/test/java', 'org.apache.commons.math3.linear', 'RealVectorFormatAbstractTest', 'testParseManyComponents', 318, 322]

    @Test
    public void testParseManyComponents() {
        RealMatrix parsed = realMatrixFormat.parse("{{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}}");
        Assert.assertEquals(24, parsed.getColumnDimension());
    }

    @Test
    public void testParseManyComponents() {
        ArrayRealVector parsed = realVectorFormat.parse("{0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0}");
        Assert.assertEquals(24, parsed.getDimension());
    }

[2548, 'src/test/java', 'org.apache.commons.math3.complex', 'ComplexTest', 'testConjugate', 159, 165]

[2548, 'src/test/java', 'org.apache.commons.math3.complex', 'ComplexTest', 'testNegate', 420, 426]

    @Test
    public void testConjugate() {
        Complex x = new Complex(3.0, 4.0);
        Complex z = x.conjugate();
        Assert.assertEquals(3.0, z.getReal(), 1.0e-5);
        Assert.assertEquals(-4.0, z.getImaginary(), 1.0e-5);
    }

    @Test
    public void testNegate() {
        Complex x = new Complex(3.0, 4.0);
        Complex z = x.negate();
        Assert.assertEquals(-3.0, z.getReal(), 1.0e-5);
        Assert.assertEquals(-4.0, z.getImaginary(), 1.0e-5);
    }

[2584, 'src/test/java', 'org.apache.commons.math3.distribution.fitting', 'MultivariateNormalMixtureExpectationMaximizationTest', 'testMaxIterationsPositive', 66, 77]

[2584, 'src/test/java', 'org.apache.commons.math3.distribution.fitting', 'MultivariateNormalMixtureExpectationMaximizationTest', 'testConvergenceException', 93, 105]

    @Test(expected = NotStrictlyPositiveException.class)
    public void testMaxIterationsPositive() {
        // Maximum iterations for fit must be positive integer
        double[][] data = getTestSamples();
        MultivariateNormalMixtureExpectationMaximization fitter =
                new MultivariateNormalMixtureExpectationMaximization(data);

        MixtureMultivariateNormalDistribution
            initialMix = MultivariateNormalMixtureExpectationMaximization.estimate(data, 2);

        fitter.fit(initialMix, 0, 1E-5);
    }

    @Test(expected = ConvergenceException.class)
    public void testConvergenceException() {
        // ConvergenceException thrown if fit terminates before threshold met
        double[][] data = getTestSamples();
        MultivariateNormalMixtureExpectationMaximization fitter
            = new MultivariateNormalMixtureExpectationMaximization(data);

        MixtureMultivariateNormalDistribution
            initialMix = MultivariateNormalMixtureExpectationMaximization.estimate(data, 2);

        // 5 iterations not enough to meet convergence threshold
        fitter.fit(initialMix, 5, 1E-5);
    }

[2614, 'src/test/java', 'org.apache.commons.math3.complex', 'ComplexTest', 'testFloatingPointEqualsPrecondition1', 499, 502]

[2614, 'src/test/java', 'org.apache.commons.math3.complex', 'ComplexTest', 'testFloatingPointEqualsPrecondition2', 503, 506]

    @Test(expected=NullPointerException.class)
    public void testFloatingPointEqualsPrecondition1() {
        Complex.equals(new Complex(3.0, 4.0), null, 3);
    }

    @Test(expected=NullPointerException.class)
    public void testFloatingPointEqualsPrecondition2() {
        Complex.equals(null, new Complex(3.0, 4.0), 3);
    }

[2620, 'src/test/java', 'org.apache.commons.math3.distribution', 'CauchyDistributionTest', 'testMedian', 80, 84]

[2620, 'src/test/java', 'org.apache.commons.math3.distribution', 'CauchyDistributionTest', 'testScale', 86, 90]

    @Test
    public void testMedian() {
        CauchyDistribution distribution = (CauchyDistribution) getDistribution();
        Assert.assertEquals(1.2, distribution.getMedian(), 0.0);
    }

    @Test
    public void testScale() {
        CauchyDistribution distribution = (CauchyDistribution) getDistribution();
        Assert.assertEquals(2.1, distribution.getScale(), 0.0);
    }

[2623, 'src/test/java', 'org.apache.commons.math3.linear', 'RealVectorAbstractTest', 'testSetEntryInvalidIndex2', 257, 260]

[2623, 'src/test/java', 'org.apache.commons.math3.linear', 'RealVectorAbstractTest', 'testAddToEntryInvalidIndex2', 304, 307]

    @Test(expected=OutOfRangeException.class)
    public void testSetEntryInvalidIndex2() {
        create(new double[4]).setEntry(4, getPreferredEntryValue());
    }

    @Test(expected=OutOfRangeException.class)
    public void testAddToEntryInvalidIndex2() {
        create(new double[3]).addToEntry(4, getPreferredEntryValue());
    }

[2636, 'src/test/java', 'org.apache.commons.math3.linear', 'Array2DRowRealMatrixTest', 'testSetRow', 813, 831]

[2636, 'src/test/java', 'org.apache.commons.math3.linear', 'BlockRealMatrixTest', 'testSetRow', 932, 950]

    @Test
    public void testSetRow() {
        RealMatrix m = new Array2DRowRealMatrix(subTestData);
        Assert.assertTrue(subRow3[0][0] != m.getRow(0)[0]);
        m.setRow(0, subRow3[0]);
        checkArrays(subRow3[0], m.getRow(0));
        try {
            m.setRow(-1, subRow3[0]);
            Assert.fail("Expecting OutOfRangeException");
        } catch (OutOfRangeException ex) {
            // expected
        }
        try {
            m.setRow(0, new double[5]);
            Assert.fail("Expecting MatrixDimensionMismatchException");
        } catch (MatrixDimensionMismatchException ex) {
            // expected
        }
    }

    @Test
    public void testSetRow() {
        RealMatrix m = new BlockRealMatrix(subTestData);
        Assert.assertTrue(subRow3[0][0] != m.getRow(0)[0]);
        m.setRow(0, subRow3[0]);
        checkArrays(subRow3[0], m.getRow(0));
        try {
            m.setRow(-1, subRow3[0]);
            Assert.fail("Expecting OutOfRangeException");
        } catch (OutOfRangeException ex) {
            // expected
        }
        try {
            m.setRow(0, new double[5]);
            Assert.fail("Expecting MatrixDimensionMismatchException");
        } catch (MatrixDimensionMismatchException ex) {
            // expected
        }
    }

[2644, 'src/test/java', 'org.apache.commons.math3.util', 'OpenIntToFieldTest', 'testPutAndGetWith0ExpectedSize', 85, 89]

[2644, 'src/test/java', 'org.apache.commons.math3.util', 'OpenIntToFieldTest', 'testPutAndGetWithExpectedSize', 91, 95]

    @Test
    public void testPutAndGetWith0ExpectedSize() {
        OpenIntToFieldHashMap<Fraction> map = new OpenIntToFieldHashMap<Fraction>(field,0);
        assertPutAndGet(map);
    }

    @Test
    public void testPutAndGetWithExpectedSize() {
        OpenIntToFieldHashMap<Fraction> map = new OpenIntToFieldHashMap<Fraction>(field,500);
        assertPutAndGet(map);
    }

[2647, 'src/test/java', 'org.apache.commons.math3.fraction', 'BigFractionFormatTest', 'setUp', 41, 45]

[2647, 'src/test/java', 'org.apache.commons.math3.fraction', 'FractionFormatTest', 'setUp', 39, 43]

    @Before
    public void setUp() {
        properFormat = BigFractionFormat.getProperInstance(getLocale());
        improperFormat = BigFractionFormat.getImproperInstance(getLocale());
    }

    @Before
    public void setUp() {
        properFormat = FractionFormat.getProperInstance(getLocale());
        improperFormat = FractionFormat.getImproperInstance(getLocale());
    }

[2648, 'src/test/java', 'org.apache.commons.math3.geometry.euclidean.threed', 'SubLineTest', 'testNoSegments', 61, 67]

[2648, 'src/test/java', 'org.apache.commons.math3.geometry.euclidean.threed', 'SubLineTest', 'testSeveralSegments', 69, 76]

    @Test
    public void testNoSegments() throws MathIllegalArgumentException {
        SubLine empty = new SubLine(new Line(new Vector3D(-1, -7, 2), new Vector3D(7, -1, 0), 1.0e-10),
                                    (IntervalsSet) new RegionFactory<Euclidean1D>().getComplement(new IntervalsSet(1.0e-10)));
        List<Segment> segments = empty.getSegments();
        Assert.assertEquals(0, segments.size());
    }

    @Test
    public void testSeveralSegments() throws MathIllegalArgumentException {
        SubLine twoSubs = new SubLine(new Line(new Vector3D(-1, -7, 2), new Vector3D(7, -1, 0), 1.0e-10),
                                      (IntervalsSet) new RegionFactory<Euclidean1D>().union(new IntervalsSet(1, 2, 1.0e-10),
                                                                                            new IntervalsSet(3, 4, 1.0e-10)));
        List<Segment> segments = twoSubs.getSegments();
        Assert.assertEquals(2, segments.size());
    }

[2655, 'src/test/java', 'org.apache.commons.math3.genetics', 'OrderedCrossoverTest', 'testCrossoverInvalidFixedLengthChromosomeFirst', 73, 86]

[2655, 'src/test/java', 'org.apache.commons.math3.genetics', 'OrderedCrossoverTest', 'testCrossoverInvalidFixedLengthChromosomeSecond', 88, 101]

    @Test(expected = MathIllegalArgumentException.class)
    public void testCrossoverInvalidFixedLengthChromosomeFirst() {
        final Integer[] p1 = new Integer[] { 1, 0, 1, 0, 0, 1, 0, 1, 1 };
        final BinaryChromosome p1c = new DummyBinaryChromosome(p1);
        final Chromosome p2c = new Chromosome() {
            public double fitness() {
                // Not important
                return 0;
            }
        };

        final CrossoverPolicy cp = new OrderedCrossover<Integer>();
        cp.crossover(p1c, p2c);
    }

    @Test(expected = MathIllegalArgumentException.class)
    public void testCrossoverInvalidFixedLengthChromosomeSecond() {
        final Integer[] p1 = new Integer[] { 1, 0, 1, 0, 0, 1, 0, 1, 1 };
        final BinaryChromosome p2c = new DummyBinaryChromosome(p1);
        final Chromosome p1c = new Chromosome() {
            public double fitness() {
                // Not important
                return 0;
            }
        };

        final CrossoverPolicy cp = new OrderedCrossover<Integer>();
        cp.crossover(p1c, p2c);
    }

[2669, 'src/test/java', 'org.apache.commons.math3.analysis.function', 'StepFunctionTest', 'testPreconditions3', 45, 48]

[2669, 'src/test/java', 'org.apache.commons.math3.analysis.function', 'StepFunctionTest', 'testPreconditions4', 50, 53]

    @Test(expected=NoDataException.class)
    public void testPreconditions3() {
        new StepFunction(new double[] {0}, new double[] {});
    }

    @Test(expected=NoDataException.class)
    public void testPreconditions4() {
        new StepFunction(new double[] {}, new double[] {0});
    }

[2671, 'src/test/java', 'org.apache.commons.math3.geometry.euclidean.twod', 'SubLineTest', 'testNoSegments', 56, 62]

[2671, 'src/test/java', 'org.apache.commons.math3.geometry.euclidean.twod', 'SubLineTest', 'testSeveralSegments', 64, 71]

    @Test
    public void testNoSegments() {
        SubLine empty = new SubLine(new Line(new Vector2D(-1, -7), new Vector2D(7, -1), 1.0e-10),
                                    new RegionFactory<Euclidean1D>().getComplement(new IntervalsSet(1.0e-10)));
        List<Segment> segments = empty.getSegments();
        Assert.assertEquals(0, segments.size());
    }

    @Test
    public void testSeveralSegments() {
        SubLine twoSubs = new SubLine(new Line(new Vector2D(-1, -7), new Vector2D(7, -1), 1.0e-10),
                                    new RegionFactory<Euclidean1D>().union(new IntervalsSet(1, 2, 1.0e-10),
                                                                           new IntervalsSet(3, 4, 1.0e-10)));
        List<Segment> segments = twoSubs.getSegments();
        Assert.assertEquals(2, segments.size());
    }

[2678, 'src/test/java', 'org.apache.commons.math3.exception', 'NumberIsTooLargeExceptionTest', 'testAccessors', 27, 33]

[2678, 'src/test/java', 'org.apache.commons.math3.exception', 'NumberIsTooSmallExceptionTest', 'testAccessors', 27, 33]

    @Test
    public void testAccessors() {
        final NumberIsTooLargeException e = new NumberIsTooLargeException(1, 0, true);
        Assert.assertEquals(1, e.getArgument());
        Assert.assertEquals(0, e.getMax());
        Assert.assertTrue(e.getBoundIsAllowed());
    }

    @Test
    public void testAccessors() {
        final NumberIsTooSmallException e = new NumberIsTooSmallException(0, 0, false);
        Assert.assertEquals(0, e.getArgument());
        Assert.assertEquals(0, e.getMin());
        Assert.assertFalse(e.getBoundIsAllowed());
    }

[2695, 'src/test/java', 'org.apache.commons.math3.complex', 'ComplexTest', 'testPowNaNBase', 886, 890]

[2695, 'src/test/java', 'org.apache.commons.math3.complex', 'ComplexTest', 'testPowNaNExponent', 892, 896]

    @Test
    public void testPowNaNBase() {
        Complex x = new Complex(3, 4);
        Assert.assertTrue(Complex.NaN.pow(x).isNaN());
    }

    @Test
    public void testPowNaNExponent() {
        Complex x = new Complex(3, 4);
        Assert.assertTrue(x.pow(Complex.NaN).isNaN());
    }

[2700, 'src/test/java', 'org.apache.commons.math3.stat', 'FrequencyTest', 'testModeDoubleNan', 413, 436]

[2700, 'src/test/java', 'org.apache.commons.math3.stat', 'FrequencyTest', 'testModeFloatNan', 438, 461]

    @Test
    public void testModeDoubleNan() {
        List<Comparable<?>> mode;
        f.addValue(Double.valueOf(Double.NaN));
        f.addValue(Double.valueOf(Double.NaN));
        f.addValue(Double.valueOf(Double.NaN));
        f.addValue(Double.valueOf(Double.NEGATIVE_INFINITY));
        f.addValue(Double.valueOf(Double.POSITIVE_INFINITY));
        f.addValue(Double.valueOf(Double.NEGATIVE_INFINITY));
        f.addValue(Double.valueOf(Double.POSITIVE_INFINITY));
        f.addValue(Double.valueOf(Double.NEGATIVE_INFINITY));
        f.addValue(Double.valueOf(Double.POSITIVE_INFINITY));
        mode = f.getMode();
        Assert.assertEquals(3, mode.size());
        Assert.assertEquals(Double.valueOf(Double.NEGATIVE_INFINITY), mode.get(0));
        Assert.assertEquals(Double.valueOf(Double.POSITIVE_INFINITY), mode.get(1));
        Assert.assertEquals(Double.valueOf(Double.NaN), mode.get(2));
        try {
            f.addValue(Float.valueOf(Float.NaN));
            Assert.fail("Expected MathIllegalArgumentException");
        } catch (MathIllegalArgumentException e) {
            // expected
        }
    }

    @Test
    public void testModeFloatNan() {
        List<Comparable<?>> mode;
        f.addValue(Float.valueOf(Float.NaN));
        f.addValue(Float.valueOf(Float.NaN));
        f.addValue(Float.valueOf(Float.NaN));
        f.addValue(Float.valueOf(Float.NEGATIVE_INFINITY));
        f.addValue(Float.valueOf(Float.POSITIVE_INFINITY));
        f.addValue(Float.valueOf(Float.NEGATIVE_INFINITY));
        f.addValue(Float.valueOf(Float.POSITIVE_INFINITY));
        f.addValue(Float.valueOf(Float.NEGATIVE_INFINITY));
        f.addValue(Float.valueOf(Float.POSITIVE_INFINITY));
        mode = f.getMode();
        Assert.assertEquals(3, mode.size());
        Assert.assertEquals(Float.valueOf(Float.NEGATIVE_INFINITY), mode.get(0));
        Assert.assertEquals(Float.valueOf(Float.POSITIVE_INFINITY), mode.get(1));
        Assert.assertEquals(Float.valueOf(Float.NaN), mode.get(2));
        try {
            f.addValue(Double.valueOf(Double.NaN));
            Assert.fail("Expected MathIllegalArgumentException");
        } catch (MathIllegalArgumentException e) {
            // expected
        }
    }

[2706, 'src/test/java', 'org.apache.commons.math3.linear', 'Array2DRowRealMatrixTest', 'testNorm', 154, 160]

[2706, 'src/test/java', 'org.apache.commons.math3.linear', 'BlockRealMatrixTest', 'testNorm', 155, 161]

    @Test
    public void testNorm() {
        Array2DRowRealMatrix m = new Array2DRowRealMatrix(testData);
        Array2DRowRealMatrix m2 = new Array2DRowRealMatrix(testData2);
        Assert.assertEquals("testData norm",14d,m.getNorm(),entryTolerance);
        Assert.assertEquals("testData2 norm",7d,m2.getNorm(),entryTolerance);
    }

    @Test
    public void testNorm() {
        BlockRealMatrix m = new BlockRealMatrix(testData);
        BlockRealMatrix m2 = new BlockRealMatrix(testData2);
        Assert.assertEquals("testData norm",14d,m.getNorm(),entryTolerance);
        Assert.assertEquals("testData2 norm",7d,m2.getNorm(),entryTolerance);
    }

[2718, 'src/test/java', 'org.apache.commons.math3.linear', 'BlockRealMatrixTest', 'testOperate', 326, 338]

[2718, 'src/test/java', 'org.apache.commons.math3.linear', 'BlockRealMatrixTest', 'testPremultiplyVector', 394, 407]

    @Test
    public void testOperate() {
        RealMatrix m = new BlockRealMatrix(id);
        assertClose(testVector, m.operate(testVector), entryTolerance);
        assertClose(testVector, m.operate(new ArrayRealVector(testVector)).toArray(), entryTolerance);
        m = new BlockRealMatrix(bigSingular);
        try {
            m.operate(testVector);
            Assert.fail("Expecting illegalArgumentException");
        } catch (MathIllegalArgumentException ex) {
            // ignored
        }
    }

    @Test
    public void testPremultiplyVector() {
        RealMatrix m = new BlockRealMatrix(testData);
        assertClose(m.preMultiply(testVector), preMultTest, normTolerance);
        assertClose(m.preMultiply(new ArrayRealVector(testVector).toArray()),
                    preMultTest, normTolerance);
        m = new BlockRealMatrix(bigSingular);
        try {
            m.preMultiply(testVector);
            Assert.fail("expecting MathIllegalArgumentException");
        } catch (MathIllegalArgumentException ex) {
            // ignored
        }
    }

[2743, 'src/test/java', 'org.apache.commons.math3.stat', 'StatUtilsTest', 'testMax', 302, 340]

[2743, 'src/test/java', 'org.apache.commons.math3.stat', 'StatUtilsTest', 'testMin', 342, 380]

    @Test
    public void testMax() {
        double[] x = null;

        try {
            StatUtils.max(x, 0, 4);
            Assert.fail("null is not a valid data array.");
        } catch (MathIllegalArgumentException ex) {
            // success
        }

        // test empty
        x = new double[] {};
        TestUtils.assertEquals(Double.NaN, StatUtils.max(x, 0, 0), TOLERANCE);

        // test one
        x = new double[] {TWO};
        TestUtils.assertEquals(TWO, StatUtils.max(x, 0, 1), TOLERANCE);

        // test many
        x = new double[] {ONE, TWO, TWO, THREE};
        TestUtils.assertEquals(THREE, StatUtils.max(x, 1, 3), TOLERANCE);

        // test first nan is ignored
        x = new double[] {NAN, TWO, THREE};
        TestUtils.assertEquals(THREE, StatUtils.max(x), TOLERANCE);

        // test middle nan is ignored
        x = new double[] {ONE, NAN, THREE};
        TestUtils.assertEquals(THREE, StatUtils.max(x), TOLERANCE);

        // test last nan is ignored
        x = new double[] {ONE, TWO, NAN};
        TestUtils.assertEquals(TWO, StatUtils.max(x), TOLERANCE);

        // test all nan returns nan
        x = new double[] {NAN, NAN, NAN};
        TestUtils.assertEquals(NAN, StatUtils.max(x), TOLERANCE);
    }

    @Test
    public void testMin() {
        double[] x = null;

        try {
            StatUtils.min(x, 0, 4);
            Assert.fail("null is not a valid data array.");
        } catch (MathIllegalArgumentException ex) {
            // success
        }

        // test empty
        x = new double[] {};
        TestUtils.assertEquals(Double.NaN, StatUtils.min(x, 0, 0), TOLERANCE);

        // test one
        x = new double[] {TWO};
        TestUtils.assertEquals(TWO, StatUtils.min(x, 0, 1), TOLERANCE);

        // test many
        x = new double[] {ONE, TWO, TWO, THREE};
        TestUtils.assertEquals(TWO, StatUtils.min(x, 1, 3), TOLERANCE);

        // test first nan is ignored
        x = new double[] {NAN, TWO, THREE};
        TestUtils.assertEquals(TWO, StatUtils.min(x), TOLERANCE);

        // test middle nan is ignored
        x = new double[] {ONE, NAN, THREE};
        TestUtils.assertEquals(ONE, StatUtils.min(x), TOLERANCE);

        // test last nan is ignored
        x = new double[] {ONE, TWO, NAN};
        TestUtils.assertEquals(ONE, StatUtils.min(x), TOLERANCE);

        // test all nan returns nan
        x = new double[] {NAN, NAN, NAN};
        TestUtils.assertEquals(NAN, StatUtils.min(x), TOLERANCE);
    }

[2753, 'src/test/java', 'org.apache.commons.math3.fraction', 'BigFractionTest', 'testGetReducedFraction', 585, 598]

[2753, 'src/test/java', 'org.apache.commons.math3.fraction', 'FractionTest', 'testGetReducedFraction', 593, 608]

    @Test
    public void testGetReducedFraction() {
        BigFraction threeFourths = new BigFraction(3, 4);
        Assert.assertTrue(threeFourths.equals(BigFraction.getReducedFraction(6, 8)));
        Assert.assertTrue(BigFraction.ZERO.equals(BigFraction.getReducedFraction(0, -1)));
        try {
            BigFraction.getReducedFraction(1, 0);
            Assert.fail("expecting ZeroException");
        } catch (ZeroException ex) {
            // expected
        }
        Assert.assertEquals(BigFraction.getReducedFraction(2, Integer.MIN_VALUE).getNumeratorAsInt(), -1);
        Assert.assertEquals(BigFraction.getReducedFraction(1, -1).getNumeratorAsInt(), -1);
    }

    @Test
    public void testGetReducedFraction() {
        Fraction threeFourths = new Fraction(3, 4);
        Assert.assertTrue(threeFourths.equals(Fraction.getReducedFraction(6, 8)));
        Assert.assertTrue(Fraction.ZERO.equals(Fraction.getReducedFraction(0, -1)));
        try {
            Fraction.getReducedFraction(1, 0);
            Assert.fail("expecting MathArithmeticException");
        } catch (MathArithmeticException ex) {
            // expected
        }
        Assert.assertEquals(Fraction.getReducedFraction
                (2, Integer.MIN_VALUE).getNumerator(),-1);
        Assert.assertEquals(Fraction.getReducedFraction
                (1, -1).getNumerator(), -1);
    }

[2782, 'src/test/java', 'org.apache.commons.math3.analysis.interpolation', 'BicubicInterpolatingFunctionTest', 'testPreconditions', 38, 113]

[2782, 'src/test/java', 'org.apache.commons.math3.analysis.interpolation', 'BicubicSplineInterpolatingFunctionTest', 'testPreconditions', 41, 116]

    @Test
    public void testPreconditions() {
        double[] xval = new double[] {3, 4, 5, 6.5};
        double[] yval = new double[] {-4, -3, -1, 2.5};
        double[][] zval = new double[xval.length][yval.length];

        @SuppressWarnings("unused")
        BivariateFunction bcf = new BicubicInterpolatingFunction(xval, yval, zval,
                                                                 zval, zval, zval);

        double[] wxval = new double[] {3, 2, 5, 6.5};
        try {
            bcf = new BicubicInterpolatingFunction(wxval, yval, zval, zval, zval, zval);
            Assert.fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException e) {
            // Expected
        }
        double[] wyval = new double[] {-4, -1, -1, 2.5};
        try {
            bcf = new BicubicInterpolatingFunction(xval, wyval, zval, zval, zval, zval);
            Assert.fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException e) {
            // Expected
        }
        double[][] wzval = new double[xval.length][yval.length - 1];
        try {
            bcf = new BicubicInterpolatingFunction(xval, yval, wzval, zval, zval, zval);
            Assert.fail("an exception should have been thrown");
        } catch (DimensionMismatchException e) {
            // Expected
        }
        try {
            bcf = new BicubicInterpolatingFunction(xval, yval, zval, wzval, zval, zval);
            Assert.fail("an exception should have been thrown");
        } catch (DimensionMismatchException e) {
            // Expected
        }
        try {
            bcf = new BicubicInterpolatingFunction(xval, yval, zval, zval, wzval, zval);
            Assert.fail("an exception should have been thrown");
        } catch (DimensionMismatchException e) {
            // Expected
        }
        try {
            bcf = new BicubicInterpolatingFunction(xval, yval, zval, zval, zval, wzval);
            Assert.fail("an exception should have been thrown");
        } catch (DimensionMismatchException e) {
            // Expected
        }

        wzval = new double[xval.length - 1][yval.length];
        try {
            bcf = new BicubicInterpolatingFunction(xval, yval, wzval, zval, zval, zval);
            Assert.fail("an exception should have been thrown");
        } catch (DimensionMismatchException e) {
            // Expected
        }
        try {
            bcf = new BicubicInterpolatingFunction(xval, yval, zval, wzval, zval, zval);
            Assert.fail("an exception should have been thrown");
        } catch (DimensionMismatchException e) {
            // Expected
        }
        try {
            bcf = new BicubicInterpolatingFunction(xval, yval, zval, zval, wzval, zval);
            Assert.fail("an exception should have been thrown");
        } catch (DimensionMismatchException e) {
            // Expected
        }
        try {
            bcf = new BicubicInterpolatingFunction(xval, yval, zval, zval, zval, wzval);
            Assert.fail("an exception should have been thrown");
        } catch (DimensionMismatchException e) {
            // Expected
        }
    }

    @Test
    public void testPreconditions() {
        double[] xval = new double[] {3, 4, 5, 6.5};
        double[] yval = new double[] {-4, -3, -1, 2.5};
        double[][] zval = new double[xval.length][yval.length];

        @SuppressWarnings("unused")
        BivariateFunction bcf = new BicubicSplineInterpolatingFunction(xval, yval, zval,
                                                                           zval, zval, zval);

        double[] wxval = new double[] {3, 2, 5, 6.5};
        try {
            bcf = new BicubicSplineInterpolatingFunction(wxval, yval, zval, zval, zval, zval);
            Assert.fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException e) {
            // Expected
        }
        double[] wyval = new double[] {-4, -1, -1, 2.5};
        try {
            bcf = new BicubicSplineInterpolatingFunction(xval, wyval, zval, zval, zval, zval);
            Assert.fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException e) {
            // Expected
        }
        double[][] wzval = new double[xval.length][yval.length - 1];
        try {
            bcf = new BicubicSplineInterpolatingFunction(xval, yval, wzval, zval, zval, zval);
            Assert.fail("an exception should have been thrown");
        } catch (DimensionMismatchException e) {
            // Expected
        }
        try {
            bcf = new BicubicSplineInterpolatingFunction(xval, yval, zval, wzval, zval, zval);
            Assert.fail("an exception should have been thrown");
        } catch (DimensionMismatchException e) {
            // Expected
        }
        try {
            bcf = new BicubicSplineInterpolatingFunction(xval, yval, zval, zval, wzval, zval);
            Assert.fail("an exception should have been thrown");
        } catch (DimensionMismatchException e) {
            // Expected
        }
        try {
            bcf = new BicubicSplineInterpolatingFunction(xval, yval, zval, zval, zval, wzval);
            Assert.fail("an exception should have been thrown");
        } catch (DimensionMismatchException e) {
            // Expected
        }

        wzval = new double[xval.length - 1][yval.length];
        try {
            bcf = new BicubicSplineInterpolatingFunction(xval, yval, wzval, zval, zval, zval);
            Assert.fail("an exception should have been thrown");
        } catch (DimensionMismatchException e) {
            // Expected
        }
        try {
            bcf = new BicubicSplineInterpolatingFunction(xval, yval, zval, wzval, zval, zval);
            Assert.fail("an exception should have been thrown");
        } catch (DimensionMismatchException e) {
            // Expected
        }
        try {
            bcf = new BicubicSplineInterpolatingFunction(xval, yval, zval, zval, wzval, zval);
            Assert.fail("an exception should have been thrown");
        } catch (DimensionMismatchException e) {
            // Expected
        }
        try {
            bcf = new BicubicSplineInterpolatingFunction(xval, yval, zval, zval, zval, wzval);
            Assert.fail("an exception should have been thrown");
        } catch (DimensionMismatchException e) {
            // Expected
        }
    }

[2786, 'src/test/java', 'org.apache.commons.math3.optimization.linear', 'SimplexSolverTest', 'testMath272', 181, 196]

[2786, 'src/test/java', 'org.apache.commons.math3.optimization.linear', 'SimplexSolverTest', 'testEpsilon', 422, 437]

    @Test
    public void testMath272() {
        LinearObjectiveFunction f = new LinearObjectiveFunction(new double[] { 2, 2, 1 }, 0);
        Collection<LinearConstraint> constraints = new ArrayList<LinearConstraint>();
        constraints.add(new LinearConstraint(new double[] { 1, 1, 0 }, Relationship.GEQ,  1));
        constraints.add(new LinearConstraint(new double[] { 1, 0, 1 }, Relationship.GEQ,  1));
        constraints.add(new LinearConstraint(new double[] { 0, 1, 0 }, Relationship.GEQ,  1));

        SimplexSolver solver = new SimplexSolver();
        PointValuePair solution = solver.optimize(f, constraints, GoalType.MINIMIZE, true);

        Assert.assertEquals(0.0, solution.getPoint()[0], .0000001);
        Assert.assertEquals(1.0, solution.getPoint()[1], .0000001);
        Assert.assertEquals(1.0, solution.getPoint()[2], .0000001);
        Assert.assertEquals(3.0, solution.getValue(), .0000001);
    }

    @Test
    public void testEpsilon() {
      LinearObjectiveFunction f =
          new LinearObjectiveFunction(new double[] { 10, 5, 1 }, 0);
      Collection<LinearConstraint> constraints = new ArrayList<LinearConstraint>();
      constraints.add(new LinearConstraint(new double[] {  9, 8, 0 }, Relationship.EQ,  17));
      constraints.add(new LinearConstraint(new double[] {  0, 7, 8 }, Relationship.LEQ,  7));
      constraints.add(new LinearConstraint(new double[] { 10, 0, 2 }, Relationship.LEQ, 10));

      SimplexSolver solver = new SimplexSolver();
      PointValuePair solution = solver.optimize(f, constraints, GoalType.MAXIMIZE, false);
      Assert.assertEquals(1.0, solution.getPoint()[0], 0.0);
      Assert.assertEquals(1.0, solution.getPoint()[1], 0.0);
      Assert.assertEquals(0.0, solution.getPoint()[2], 0.0);
      Assert.assertEquals(15.0, solution.getValue(), 0.0);
  }

[2794, 'src/test/java', 'org.apache.commons.math3.stat', 'StatUtilsTest', 'testVariance', 245, 271]

[2794, 'src/test/java', 'org.apache.commons.math3.stat', 'StatUtilsTest', 'testPopulationVariance', 273, 299]

    @Test
    public void testVariance() {
        double[] x = null;

        try {
            StatUtils.variance(x, 0, 4);
            Assert.fail("null is not a valid data array.");
        } catch (MathIllegalArgumentException ex) {
            // success
        }

        // test empty
        x = new double[] {};
        TestUtils.assertEquals(Double.NaN, StatUtils.variance(x, 0, 0), TOLERANCE);

        // test one
        x = new double[] {TWO};
        TestUtils.assertEquals(0.0, StatUtils.variance(x, 0, 1), TOLERANCE);

        // test many
        x = new double[] {ONE, TWO, TWO, THREE};
        TestUtils.assertEquals(0.5, StatUtils.variance(x, 2, 2), TOLERANCE);

        // test precomputed mean
        x = new double[] {ONE, TWO, TWO, THREE};
        TestUtils.assertEquals(0.5, StatUtils.variance(x,2.5, 2, 2), TOLERANCE);
    }

    @Test
    public void testPopulationVariance() {
        double[] x = null;

        try {
            StatUtils.variance(x, 0, 4);
            Assert.fail("null is not a valid data array.");
        } catch (MathIllegalArgumentException ex) {
            // success
        }

        // test empty
        x = new double[] {};
        TestUtils.assertEquals(Double.NaN, StatUtils.populationVariance(x, 0, 0), TOLERANCE);

        // test one
        x = new double[] {TWO};
        TestUtils.assertEquals(0.0, StatUtils.populationVariance(x, 0, 1), TOLERANCE);

        // test many
        x = new double[] {ONE, TWO, TWO, THREE};
        TestUtils.assertEquals(0.25, StatUtils.populationVariance(x, 0, 2), TOLERANCE);

        // test precomputed mean
        x = new double[] {ONE, TWO, TWO, THREE};
        TestUtils.assertEquals(0.25, StatUtils.populationVariance(x, 2.5, 2, 2), TOLERANCE);
    }

[2801, 'src/test/java', 'org.apache.commons.math3.optim.linear', 'SimplexSolverTest', 'testModelWithNoArtificialVars', 501, 515]

[2801, 'src/test/java', 'org.apache.commons.math3.optim.linear', 'SimplexSolverTest', 'testSimplexSolver', 467, 482]

    @Test
    public void testModelWithNoArtificialVars() {
        LinearObjectiveFunction f = new LinearObjectiveFunction(new double[] { 15, 10 }, 0);
        Collection<LinearConstraint> constraints = new ArrayList<LinearConstraint>();
        constraints.add(new LinearConstraint(new double[] { 1, 0 }, Relationship.LEQ, 2));
        constraints.add(new LinearConstraint(new double[] { 0, 1 }, Relationship.LEQ, 3));
        constraints.add(new LinearConstraint(new double[] { 1, 1 }, Relationship.LEQ, 4));

        SimplexSolver solver = new SimplexSolver();
        PointValuePair solution = solver.optimize(DEFAULT_MAX_ITER, f, new LinearConstraintSet(constraints),
                                                  GoalType.MAXIMIZE, new NonNegativeConstraint(false));
        Assert.assertEquals(2.0, solution.getPoint()[0], 0.0);
        Assert.assertEquals(2.0, solution.getPoint()[1], 0.0);
        Assert.assertEquals(50.0, solution.getValue(), 0.0);
    }

    @Test
    public void testSimplexSolver() {
        LinearObjectiveFunction f =
            new LinearObjectiveFunction(new double[] { 15, 10 }, 7);
        Collection<LinearConstraint> constraints = new ArrayList<LinearConstraint>();
        constraints.add(new LinearConstraint(new double[] { 1, 0 }, Relationship.LEQ, 2));
        constraints.add(new LinearConstraint(new double[] { 0, 1 }, Relationship.LEQ, 3));
        constraints.add(new LinearConstraint(new double[] { 1, 1 }, Relationship.EQ, 4));

        SimplexSolver solver = new SimplexSolver();
        PointValuePair solution = solver.optimize(DEFAULT_MAX_ITER, f, new LinearConstraintSet(constraints),
                                                  GoalType.MAXIMIZE, new NonNegativeConstraint(true));
        Assert.assertEquals(2.0, solution.getPoint()[0], 0.0);
        Assert.assertEquals(2.0, solution.getPoint()[1], 0.0);
        Assert.assertEquals(57.0, solution.getValue(), 0.0);
    }

[2841, 'src/test/java', 'org.apache.commons.math3.util', 'CombinationsTest', 'testLexicographicComparator', 91, 102]

[2841, 'src/test/java', 'org.apache.commons.math3.util', 'CombinationsTest', 'testLexicographicComparatorUnsorted', 107, 118]

    @Test
    public void testLexicographicComparator() {
        final int n = 5;
        final int k = 3;
        final Comparator<int[]> comp = new Combinations(n, k).comparator();
        Assert.assertEquals(1, comp.compare(new int[] {1, 2, 4},
                                            new int[] {1, 2, 3}));
        Assert.assertEquals(-1, comp.compare(new int[] {0, 1, 4},
                                             new int[] {0, 2, 4}));
        Assert.assertEquals(0, comp.compare(new int[] {1, 3, 4},
                                            new int[] {1, 3, 4}));
    }

    @Test
    public void testLexicographicComparatorUnsorted() {
        final int n = 5;
        final int k = 3;
        final Comparator<int[]> comp = new Combinations(n, k).comparator();
        Assert.assertEquals(1, comp.compare(new int[] {1, 4, 2},
                                            new int[] {1, 3, 2}));
        Assert.assertEquals(-1, comp.compare(new int[] {0, 4, 1},
                                             new int[] {0, 4, 2}));
        Assert.assertEquals(0, comp.compare(new int[] {1, 4, 3},
                                            new int[] {1, 3, 4}));
    }

[2861, 'src/test/java', 'org.apache.commons.math3.linear', 'Array2DRowRealMatrixTest', 'testAddFail', 141, 151]

[2861, 'src/test/java', 'org.apache.commons.math3.linear', 'BlockRealMatrixTest', 'testAddFail', 142, 152]

    @Test
    public void testAddFail() {
        Array2DRowRealMatrix m = new Array2DRowRealMatrix(testData);
        Array2DRowRealMatrix m2 = new Array2DRowRealMatrix(testData2);
        try {
            m.add(m2);
            Assert.fail("MathIllegalArgumentException expected");
        } catch (MathIllegalArgumentException ex) {
            // ignored
        }
    }

    @Test
    public void testAddFail() {
        BlockRealMatrix m = new BlockRealMatrix(testData);
        BlockRealMatrix m2 = new BlockRealMatrix(testData2);
        try {
            m.add(m2);
            Assert.fail("MathIllegalArgumentException expected");
        } catch (MathIllegalArgumentException ex) {
            // ignored
        }
    }

[2886, 'src/test/java', 'org.apache.commons.math3.analysis.integration', 'IterativeLegendreGaussIntegratorTest', 'testSinFunction', 34, 52]

[2886, 'src/test/java', 'org.apache.commons.math3.analysis.integration', 'LegendreGaussIntegratorTest', 'testSinFunction', 34, 51]

    @Test
    public void testSinFunction() {
        UnivariateFunction f = new Sin();
        BaseAbstractUnivariateIntegrator integrator
            = new IterativeLegendreGaussIntegrator(5, 1.0e-14, 1.0e-10, 2, 15);
        double min, max, expected, result, tolerance;

        min = 0; max = FastMath.PI; expected = 2;
        tolerance = FastMath.max(integrator.getAbsoluteAccuracy(),
                             FastMath.abs(expected * integrator.getRelativeAccuracy()));
        result = integrator.integrate(10000, f, min, max);
        Assert.assertEquals(expected, result, tolerance);

        min = -FastMath.PI/3; max = 0; expected = -0.5;
        tolerance = FastMath.max(integrator.getAbsoluteAccuracy(),
                FastMath.abs(expected * integrator.getRelativeAccuracy()));
        result = integrator.integrate(10000, f, min, max);
        Assert.assertEquals(expected, result, tolerance);
    }

    @Test
    public void testSinFunction() {
        UnivariateFunction f = new Sin();
        BaseAbstractUnivariateIntegrator integrator = new LegendreGaussIntegrator(5, 1.0e-14, 1.0e-10, 2, 15);
        double min, max, expected, result, tolerance;

        min = 0; max = FastMath.PI; expected = 2;
        tolerance = FastMath.max(integrator.getAbsoluteAccuracy(),
                             FastMath.abs(expected * integrator.getRelativeAccuracy()));
        result = integrator.integrate(10000, f, min, max);
        Assert.assertEquals(expected, result, tolerance);

        min = -FastMath.PI/3; max = 0; expected = -0.5;
        tolerance = FastMath.max(integrator.getAbsoluteAccuracy(),
                FastMath.abs(expected * integrator.getRelativeAccuracy()));
        result = integrator.integrate(10000, f, min, max);
        Assert.assertEquals(expected, result, tolerance);
    }

[2902, 'src/test/java', 'org.apache.commons.math3.analysis.interpolation', 'TricubicInterpolatingFunctionTest', 'testPreconditions', 38, 272]

[2902, 'src/test/java', 'org.apache.commons.math3.analysis.interpolation', 'TricubicSplineInterpolatingFunctionTest', 'testPreconditions', 36, 270]

    @Test
    public void testPreconditions() {
        double[] xval = new double[] {3, 4, 5, 6.5};
        double[] yval = new double[] {-4, -3, -1, 2.5};
        double[] zval = new double[] {-12, -8, -5.5, -3, 0, 2.5};
        double[][][] fval = new double[xval.length][yval.length][zval.length];

        @SuppressWarnings("unused")
        TrivariateFunction tcf = new TricubicInterpolatingFunction(xval, yval, zval,
                                                                   fval, fval, fval, fval,
                                                                   fval, fval, fval, fval);

        double[] wxval = new double[] {3, 2, 5, 6.5};
        try {
            tcf = new TricubicInterpolatingFunction(wxval, yval, zval,
                                                    fval, fval, fval, fval,
                                                    fval, fval, fval, fval);
            Assert.fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException e) {
            // Expected
        }
        double[] wyval = new double[] {-4, -1, -1, 2.5};
        try {
            tcf = new TricubicInterpolatingFunction(xval, wyval, zval,
                                                    fval, fval, fval, fval,
                                                    fval, fval, fval, fval);
            Assert.fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException e) {
            // Expected
        }
        double[] wzval = new double[] {-12, -8, -9, -3, 0, 2.5};
        try {
            tcf = new TricubicInterpolatingFunction(xval, yval, wzval,
                                                    fval, fval, fval, fval,
                                                    fval, fval, fval, fval);
            Assert.fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException e) {
            // Expected
        }
        double[][][] wfval = new double[xval.length - 1][yval.length - 1][zval.length];
        try {
            tcf = new TricubicInterpolatingFunction(xval, yval, zval,
                                                    wfval, fval, fval, fval,
                                                    fval, fval, fval, fval);
            Assert.fail("an exception should have been thrown");
        } catch (DimensionMismatchException e) {
            // Expected
        }
        try {
            tcf = new TricubicInterpolatingFunction(xval, yval, zval,
                                                    fval, wfval, fval, fval,
                                                    fval, fval, fval, fval);
            Assert.fail("an exception should have been thrown");
        } catch (DimensionMismatchException e) {
            // Expected
        }
        try {
            tcf = new TricubicInterpolatingFunction(xval, yval, zval,
                                                    fval, fval, wfval, fval,
                                                    fval, fval, fval, fval);
            Assert.fail("an exception should have been thrown");
        } catch (DimensionMismatchException e) {
            // Expected
        }
        try {
            tcf = new TricubicInterpolatingFunction(xval, yval, zval,
                                                    fval, fval, fval, wfval,
                                                    fval, fval, fval, fval);
            Assert.fail("an exception should have been thrown");
        } catch (DimensionMismatchException e) {
            // Expected
        }
        try {
            tcf = new TricubicInterpolatingFunction(xval, yval, zval,
                                                    fval, fval, fval, fval,
                                                    wfval, fval, fval, fval);
            Assert.fail("an exception should have been thrown");
        } catch (DimensionMismatchException e) {
            // Expected
        }
        try {
            tcf = new TricubicInterpolatingFunction(xval, yval, zval,
                                                    fval, fval, fval, fval,
                                                    fval, wfval, fval, fval);
            Assert.fail("an exception should have been thrown");
        } catch (DimensionMismatchException e) {
            // Expected
        }
        try {
            tcf = new TricubicInterpolatingFunction(xval, yval, zval,
                                                    fval, fval, fval, fval,
                                                    fval, fval, wfval, fval);
            Assert.fail("an exception should have been thrown");
        } catch (DimensionMismatchException e) {
            // Expected
        }
        try {
            tcf = new TricubicInterpolatingFunction(xval, yval, zval,
                                                    fval, fval, fval, fval,
                                                    fval, fval, fval, wfval);
            Assert.fail("an exception should have been thrown");
        } catch (DimensionMismatchException e) {
            // Expected
        }
        wfval = new double[xval.length][yval.length - 1][zval.length];
        try {
            tcf = new TricubicInterpolatingFunction(xval, yval, zval,
                                                    wfval, fval, fval, fval,
                                                    fval, fval, fval, fval);
            Assert.fail("an exception should have been thrown");
        } catch (DimensionMismatchException e) {
            // Expected
        }
        try {
            tcf = new TricubicInterpolatingFunction(xval, yval, zval,
                                                    fval, wfval, fval, fval,
                                                    fval, fval, fval, fval);
            Assert.fail("an exception should have been thrown");
        } catch (DimensionMismatchException e) {
            // Expected
        }
        try {
            tcf = new TricubicInterpolatingFunction(xval, yval, zval,
                                                    fval, fval, wfval, fval,
                                                    fval, fval, fval, fval);
            Assert.fail("an exception should have been thrown");
        } catch (DimensionMismatchException e) {
            // Expected
        }
        try {
            tcf = new TricubicInterpolatingFunction(xval, yval, zval,
                                                    fval, fval, fval, wfval,
                                                    fval, fval, fval, fval);
            Assert.fail("an exception should have been thrown");
        } catch (DimensionMismatchException e) {
            // Expected
        }
        try {
            tcf = new TricubicInterpolatingFunction(xval, yval, zval,
                                                    fval, fval, fval, fval,
                                                    wfval, fval, fval, fval);
            Assert.fail("an exception should have been thrown");
        } catch (DimensionMismatchException e) {
            // Expected
        }
        try {
            tcf = new TricubicInterpolatingFunction(xval, yval, zval,
                                                    fval, fval, fval, fval,
                                                    fval, wfval, fval, fval);
            Assert.fail("an exception should have been thrown");
        } catch (DimensionMismatchException e) {
            // Expected
        }
        try {
            tcf = new TricubicInterpolatingFunction(xval, yval, zval,
                                                    fval, fval, fval, fval,
                                                    fval, fval, wfval, fval);
            Assert.fail("an exception should have been thrown");
        } catch (DimensionMismatchException e) {
            // Expected
        }
        try {
            tcf = new TricubicInterpolatingFunction(xval, yval, zval,
                                                    fval, fval, fval, fval,
                                                    fval, fval, fval, wfval);
            Assert.fail("an exception should have been thrown");
        } catch (DimensionMismatchException e) {
            // Expected
        }
        wfval = new double[xval.length][yval.length][zval.length - 1];
        try {
            tcf = new TricubicInterpolatingFunction(xval, yval, zval,
                                                    wfval, fval, fval, fval,
                                                    fval, fval, fval, fval);
            Assert.fail("an exception should have been thrown");
        } catch (DimensionMismatchException e) {
            // Expected
        }
        try {
            tcf = new TricubicInterpolatingFunction(xval, yval, zval,
                                                    fval, wfval, fval, fval,
                                                    fval, fval, fval, fval);
            Assert.fail("an exception should have been thrown");
        } catch (DimensionMismatchException e) {
            // Expected
        }
        try {
            tcf = new TricubicInterpolatingFunction(xval, yval, zval,
                                                    fval, fval, wfval, fval,
                                                    fval, fval, fval, fval);
            Assert.fail("an exception should have been thrown");
        } catch (DimensionMismatchException e) {
            // Expected
        }
        try {
            tcf = new TricubicInterpolatingFunction(xval, yval, zval,
                                                    fval, fval, fval, wfval,
                                                    fval, fval, fval, fval);
            Assert.fail("an exception should have been thrown");
        } catch (DimensionMismatchException e) {
            // Expected
        }
        try {
            tcf = new TricubicInterpolatingFunction(xval, yval, zval,
                                                    fval, fval, fval, fval,
                                                    wfval, fval, fval, fval);
            Assert.fail("an exception should have been thrown");
        } catch (DimensionMismatchException e) {
            // Expected
        }
        try {
            tcf = new TricubicInterpolatingFunction(xval, yval, zval,
                                                    fval, fval, fval, fval,
                                                    fval, wfval, fval, fval);
            Assert.fail("an exception should have been thrown");
        } catch (DimensionMismatchException e) {
            // Expected
        }
        try {
            tcf = new TricubicInterpolatingFunction(xval, yval, zval,
                                                    fval, fval, fval, fval,
                                                    fval, fval, wfval, fval);
            Assert.fail("an exception should have been thrown");
        } catch (DimensionMismatchException e) {
            // Expected
        }
        try {
            tcf = new TricubicInterpolatingFunction(xval, yval, zval,
                                                    fval, fval, fval, fval,
                                                    fval, fval, fval, wfval);
            Assert.fail("an exception should have been thrown");
        } catch (DimensionMismatchException e) {
            // Expected
        }
    }

    @Test
    public void testPreconditions() {
        double[] xval = new double[] {3, 4, 5, 6.5};
        double[] yval = new double[] {-4, -3, -1, 2.5};
        double[] zval = new double[] {-12, -8, -5.5, -3, 0, 2.5};
        double[][][] fval = new double[xval.length][yval.length][zval.length];

        @SuppressWarnings("unused")
        TrivariateFunction tcf = new TricubicSplineInterpolatingFunction(xval, yval, zval,
                                                                             fval, fval, fval, fval,
                                                                             fval, fval, fval, fval);

        double[] wxval = new double[] {3, 2, 5, 6.5};
        try {
            tcf = new TricubicSplineInterpolatingFunction(wxval, yval, zval,
                                                          fval, fval, fval, fval,
                                                          fval, fval, fval, fval);
            Assert.fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException e) {
            // Expected
        }
        double[] wyval = new double[] {-4, -1, -1, 2.5};
        try {
            tcf = new TricubicSplineInterpolatingFunction(xval, wyval, zval,
                                                          fval, fval, fval, fval,
                                                          fval, fval, fval, fval);
            Assert.fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException e) {
            // Expected
        }
        double[] wzval = new double[] {-12, -8, -9, -3, 0, 2.5};
        try {
            tcf = new TricubicSplineInterpolatingFunction(xval, yval, wzval,
                                                          fval, fval, fval, fval,
                                                          fval, fval, fval, fval);
            Assert.fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException e) {
            // Expected
        }
        double[][][] wfval = new double[xval.length - 1][yval.length - 1][zval.length];
        try {
            tcf = new TricubicSplineInterpolatingFunction(xval, yval, zval,
                                                          wfval, fval, fval, fval,
                                                          fval, fval, fval, fval);
            Assert.fail("an exception should have been thrown");
        } catch (DimensionMismatchException e) {
            // Expected
        }
        try {
            tcf = new TricubicSplineInterpolatingFunction(xval, yval, zval,
                                                          fval, wfval, fval, fval,
                                                          fval, fval, fval, fval);
            Assert.fail("an exception should have been thrown");
        } catch (DimensionMismatchException e) {
            // Expected
        }
        try {
            tcf = new TricubicSplineInterpolatingFunction(xval, yval, zval,
                                                          fval, fval, wfval, fval,
                                                          fval, fval, fval, fval);
            Assert.fail("an exception should have been thrown");
        } catch (DimensionMismatchException e) {
            // Expected
        }
        try {
            tcf = new TricubicSplineInterpolatingFunction(xval, yval, zval,
                                                          fval, fval, fval, wfval,
                                                          fval, fval, fval, fval);
            Assert.fail("an exception should have been thrown");
        } catch (DimensionMismatchException e) {
            // Expected
        }
        try {
            tcf = new TricubicSplineInterpolatingFunction(xval, yval, zval,
                                                          fval, fval, fval, fval,
                                                          wfval, fval, fval, fval);
            Assert.fail("an exception should have been thrown");
        } catch (DimensionMismatchException e) {
            // Expected
        }
        try {
            tcf = new TricubicSplineInterpolatingFunction(xval, yval, zval,
                                                          fval, fval, fval, fval,
                                                          fval, wfval, fval, fval);
            Assert.fail("an exception should have been thrown");
        } catch (DimensionMismatchException e) {
            // Expected
        }
        try {
            tcf = new TricubicSplineInterpolatingFunction(xval, yval, zval,
                                                          fval, fval, fval, fval,
                                                          fval, fval, wfval, fval);
            Assert.fail("an exception should have been thrown");
        } catch (DimensionMismatchException e) {
            // Expected
        }
        try {
            tcf = new TricubicSplineInterpolatingFunction(xval, yval, zval,
                                                          fval, fval, fval, fval,
                                                          fval, fval, fval, wfval);
            Assert.fail("an exception should have been thrown");
        } catch (DimensionMismatchException e) {
            // Expected
        }
        wfval = new double[xval.length][yval.length - 1][zval.length];
        try {
            tcf = new TricubicSplineInterpolatingFunction(xval, yval, zval,
                                                          wfval, fval, fval, fval,
                                                          fval, fval, fval, fval);
            Assert.fail("an exception should have been thrown");
        } catch (DimensionMismatchException e) {
            // Expected
        }
        try {
            tcf = new TricubicSplineInterpolatingFunction(xval, yval, zval,
                                                          fval, wfval, fval, fval,
                                                          fval, fval, fval, fval);
            Assert.fail("an exception should have been thrown");
        } catch (DimensionMismatchException e) {
            // Expected
        }
        try {
            tcf = new TricubicSplineInterpolatingFunction(xval, yval, zval,
                                                          fval, fval, wfval, fval,
                                                          fval, fval, fval, fval);
            Assert.fail("an exception should have been thrown");
        } catch (DimensionMismatchException e) {
            // Expected
        }
        try {
            tcf = new TricubicSplineInterpolatingFunction(xval, yval, zval,
                                                          fval, fval, fval, wfval,
                                                          fval, fval, fval, fval);
            Assert.fail("an exception should have been thrown");
        } catch (DimensionMismatchException e) {
            // Expected
        }
        try {
            tcf = new TricubicSplineInterpolatingFunction(xval, yval, zval,
                                                          fval, fval, fval, fval,
                                                          wfval, fval, fval, fval);
            Assert.fail("an exception should have been thrown");
        } catch (DimensionMismatchException e) {
            // Expected
        }
        try {
            tcf = new TricubicSplineInterpolatingFunction(xval, yval, zval,
                                                          fval, fval, fval, fval,
                                                          fval, wfval, fval, fval);
            Assert.fail("an exception should have been thrown");
        } catch (DimensionMismatchException e) {
            // Expected
        }
        try {
            tcf = new TricubicSplineInterpolatingFunction(xval, yval, zval,
                                                          fval, fval, fval, fval,
                                                          fval, fval, wfval, fval);
            Assert.fail("an exception should have been thrown");
        } catch (DimensionMismatchException e) {
            // Expected
        }
        try {
            tcf = new TricubicSplineInterpolatingFunction(xval, yval, zval,
                                                          fval, fval, fval, fval,
                                                          fval, fval, fval, wfval);
            Assert.fail("an exception should have been thrown");
        } catch (DimensionMismatchException e) {
            // Expected
        }
        wfval = new double[xval.length][yval.length][zval.length - 1];
        try {
            tcf = new TricubicSplineInterpolatingFunction(xval, yval, zval,
                                                          wfval, fval, fval, fval,
                                                          fval, fval, fval, fval);
            Assert.fail("an exception should have been thrown");
        } catch (DimensionMismatchException e) {
            // Expected
        }
        try {
            tcf = new TricubicSplineInterpolatingFunction(xval, yval, zval,
                                                          fval, wfval, fval, fval,
                                                          fval, fval, fval, fval);
            Assert.fail("an exception should have been thrown");
        } catch (DimensionMismatchException e) {
            // Expected
        }
        try {
            tcf = new TricubicSplineInterpolatingFunction(xval, yval, zval,
                                                          fval, fval, wfval, fval,
                                                          fval, fval, fval, fval);
            Assert.fail("an exception should have been thrown");
        } catch (DimensionMismatchException e) {
            // Expected
        }
        try {
            tcf = new TricubicSplineInterpolatingFunction(xval, yval, zval,
                                                          fval, fval, fval, wfval,
                                                          fval, fval, fval, fval);
            Assert.fail("an exception should have been thrown");
        } catch (DimensionMismatchException e) {
            // Expected
        }
        try {
            tcf = new TricubicSplineInterpolatingFunction(xval, yval, zval,
                                                          fval, fval, fval, fval,
                                                          wfval, fval, fval, fval);
            Assert.fail("an exception should have been thrown");
        } catch (DimensionMismatchException e) {
            // Expected
        }
        try {
            tcf = new TricubicSplineInterpolatingFunction(xval, yval, zval,
                                                          fval, fval, fval, fval,
                                                          fval, wfval, fval, fval);
            Assert.fail("an exception should have been thrown");
        } catch (DimensionMismatchException e) {
            // Expected
        }
        try {
            tcf = new TricubicSplineInterpolatingFunction(xval, yval, zval,
                                                          fval, fval, fval, fval,
                                                          fval, fval, wfval, fval);
            Assert.fail("an exception should have been thrown");
        } catch (DimensionMismatchException e) {
            // Expected
        }
        try {
            tcf = new TricubicSplineInterpolatingFunction(xval, yval, zval,
                                                          fval, fval, fval, fval,
                                                          fval, fval, fval, wfval);
            Assert.fail("an exception should have been thrown");
        } catch (DimensionMismatchException e) {
            // Expected
        }
    }

[2913, 'src/test/java', 'org.apache.commons.math3.linear', 'ConjugateGradientTest', 'testPreconditionedSolution', 306, 328]

[2913, 'src/test/java', 'org.apache.commons.math3.linear', 'SymmLQTest', 'testPreconditionedSolution', 407, 429]

    @Test
    public void testPreconditionedSolution() {
        final int n = 8;
        final int maxIterations = 100;
        final RealLinearOperator a = new HilbertMatrix(n);
        final InverseHilbertMatrix ainv = new InverseHilbertMatrix(n);
        final RealLinearOperator m = JacobiPreconditioner.create(a);
        final PreconditionedIterativeLinearSolver solver;
        solver = new ConjugateGradient(maxIterations, 1E-15, true);
        final RealVector b = new ArrayRealVector(n);
        for (int j = 0; j < n; j++) {
            b.set(0.);
            b.setEntry(j, 1.);
            final RealVector x = solver.solve(a, m, b);
            for (int i = 0; i < n; i++) {
                final double actual = x.getEntry(i);
                final double expected = ainv.getEntry(i, j);
                final double delta = 1E-6 * FastMath.abs(expected);
                final String msg = String.format("coefficient (%d, %d)", i, j);
                Assert.assertEquals(msg, expected, actual, delta);
            }
        }
    }

    @Test
    public void testPreconditionedSolution() {
        final int n = 8;
        final int maxIterations = 100;
        final RealLinearOperator a = new HilbertMatrix(n);
        final InverseHilbertMatrix ainv = new InverseHilbertMatrix(n);
        final RealLinearOperator m = JacobiPreconditioner.create(a);
        final PreconditionedIterativeLinearSolver solver;
        solver = new SymmLQ(maxIterations, 1E-15, true);
        final RealVector b = new ArrayRealVector(n);
        for (int j = 0; j < n; j++) {
            b.set(0.);
            b.setEntry(j, 1.);
            final RealVector x = solver.solve(a, m, b);
            for (int i = 0; i < n; i++) {
                final double actual = x.getEntry(i);
                final double expected = ainv.getEntry(i, j);
                final double delta = 1E-6 * FastMath.abs(expected);
                final String msg = String.format("coefficient (%d, %d)", i, j);
                Assert.assertEquals(msg, expected, actual, delta);
            }
        }
    }

[2935, 'src/test/java', 'org.apache.commons.math3.linear', 'SingularValueDecompositionTest', 'testMoreRows', 46, 59]

[2935, 'src/test/java', 'org.apache.commons.math3.linear', 'SingularValueDecompositionTest', 'testMoreColumns', 61, 74]

    @Test
    public void testMoreRows() {
        final double[] singularValues = { 123.456, 2.3, 1.001, 0.999 };
        final int rows    = singularValues.length + 2;
        final int columns = singularValues.length;
        Random r = new Random(15338437322523l);
        SingularValueDecomposition svd =
            new SingularValueDecomposition(createTestMatrix(r, rows, columns, singularValues));
        double[] computedSV = svd.getSingularValues();
        Assert.assertEquals(singularValues.length, computedSV.length);
        for (int i = 0; i < singularValues.length; ++i) {
            Assert.assertEquals(singularValues[i], computedSV[i], 1.0e-10);
        }
    }

    @Test
    public void testMoreColumns() {
        final double[] singularValues = { 123.456, 2.3, 1.001, 0.999 };
        final int rows    = singularValues.length;
        final int columns = singularValues.length + 2;
        Random r = new Random(732763225836210l);
        SingularValueDecomposition svd =
            new SingularValueDecomposition(createTestMatrix(r, rows, columns, singularValues));
        double[] computedSV = svd.getSingularValues();
        Assert.assertEquals(singularValues.length, computedSV.length);
        for (int i = 0; i < singularValues.length; ++i) {
            Assert.assertEquals(singularValues[i], computedSV[i], 1.0e-10);
        }
    }

[2943, 'src/test/java', 'org.apache.commons.math3.analysis.integration', 'IterativeLegendreGaussIntegratorTest', 'testQuinticFunction', 54, 76]

[2943, 'src/test/java', 'org.apache.commons.math3.analysis.integration', 'LegendreGaussIntegratorTest', 'testQuinticFunction', 53, 75]

    @Test
    public void testQuinticFunction() {
        UnivariateFunction f = new QuinticFunction();
        UnivariateIntegrator integrator =
                new IterativeLegendreGaussIntegrator(3,
                                                     BaseAbstractUnivariateIntegrator.DEFAULT_RELATIVE_ACCURACY,
                                                     BaseAbstractUnivariateIntegrator.DEFAULT_ABSOLUTE_ACCURACY,
                                                     BaseAbstractUnivariateIntegrator.DEFAULT_MIN_ITERATIONS_COUNT,
                                                     64);
        double min, max, expected, result;

        min = 0; max = 1; expected = -1.0/48;
        result = integrator.integrate(10000, f, min, max);
        Assert.assertEquals(expected, result, 1.0e-16);

        min = 0; max = 0.5; expected = 11.0/768;
        result = integrator.integrate(10000, f, min, max);
        Assert.assertEquals(expected, result, 1.0e-16);

        min = -1; max = 4; expected = 2048/3.0 - 78 + 1.0/48;
        result = integrator.integrate(10000, f, min, max);
        Assert.assertEquals(expected, result, 1.0e-16);
    }

    @Test
    public void testQuinticFunction() {
        UnivariateFunction f = new QuinticFunction();
        UnivariateIntegrator integrator =
                new LegendreGaussIntegrator(3,
                                            BaseAbstractUnivariateIntegrator.DEFAULT_RELATIVE_ACCURACY,
                                            BaseAbstractUnivariateIntegrator.DEFAULT_ABSOLUTE_ACCURACY,
                                            BaseAbstractUnivariateIntegrator.DEFAULT_MIN_ITERATIONS_COUNT,
                                            64);
        double min, max, expected, result;

        min = 0; max = 1; expected = -1.0/48;
        result = integrator.integrate(10000, f, min, max);
        Assert.assertEquals(expected, result, 1.0e-16);

        min = 0; max = 0.5; expected = 11.0/768;
        result = integrator.integrate(10000, f, min, max);
        Assert.assertEquals(expected, result, 1.0e-16);

        min = -1; max = 4; expected = 2048/3.0 - 78 + 1.0/48;
        result = integrator.integrate(10000, f, min, max);
        Assert.assertEquals(expected, result, 1.0e-16);
    }

[2960, 'src/test/java', 'org.apache.commons.math3.linear', 'BlockRealMatrixTest', 'testGetSetRowVectorLarge', 825, 842]

[2960, 'src/test/java', 'org.apache.commons.math3.linear', 'BlockRealMatrixTest', 'testGetSetColumnVectorLarge', 886, 903]

    @Test
    public void testGetSetRowVectorLarge() {
        int n = 3 * BlockRealMatrix.BLOCK_SIZE;
        RealMatrix m = new BlockRealMatrix(n, n);
        RealVector sub = new ArrayRealVector(n, 1.0);

        m.setRowVector(2, sub);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                if (i != 2) {
                    Assert.assertEquals(0.0, m.getEntry(i, j), 0.0);
                } else {
                    Assert.assertEquals(1.0, m.getEntry(i, j), 0.0);
                }
            }
        }
        Assert.assertEquals(sub, m.getRowVector(2));
    }

    @Test
    public void testGetSetColumnVectorLarge() {
        int n = 3 * BlockRealMatrix.BLOCK_SIZE;
        RealMatrix m = new BlockRealMatrix(n, n);
        RealVector sub = new ArrayRealVector(n, 1.0);

        m.setColumnVector(2, sub);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                if (j != 2) {
                    Assert.assertEquals(0.0, m.getEntry(i, j), 0.0);
                } else {
                    Assert.assertEquals(1.0, m.getEntry(i, j), 0.0);
                }
            }
        }
        Assert.assertEquals(sub, m.getColumnVector(2));
    }

[2961, 'src/test/java', 'org.apache.commons.math3.linear', 'SparseFieldVectorTest', 'testWalkInDefaultOrderChangingVisitor2', 545, 593]

[2961, 'src/test/java', 'org.apache.commons.math3.linear', 'SparseFieldVectorTest', 'testWalkInOptimizedOrderChangingVisitor2', 680, 728]

    @Test
    public void testWalkInDefaultOrderChangingVisitor2() {
        final SparseFieldVector<Fraction> v = create(5);
        final FieldVectorChangingVisitor<Fraction> visitor;
        visitor = new FieldVectorChangingVisitor<Fraction>() {

            public Fraction visit(int index, Fraction value) {
                return Fraction.ZERO;
            }

            public void start(int dimension, int start, int end) {
                // Do nothing
            }

            public Fraction end() {
                return Fraction.ZERO;
            }
        };
        try {
            v.walkInDefaultOrder(visitor, -1, 4);
            Assert.fail();
        } catch (OutOfRangeException e) {
            // Expected behavior
        }
        try {
            v.walkInDefaultOrder(visitor, 5, 4);
            Assert.fail();
        } catch (OutOfRangeException e) {
            // Expected behavior
        }
        try {
            v.walkInDefaultOrder(visitor, 0, -1);
            Assert.fail();
        } catch (OutOfRangeException e) {
            // Expected behavior
        }
        try {
            v.walkInDefaultOrder(visitor, 0, 5);
            Assert.fail();
        } catch (OutOfRangeException e) {
            // Expected behavior
        }
        try {
            v.walkInDefaultOrder(visitor, 4, 0);
            Assert.fail();
        } catch (NumberIsTooSmallException e) {
            // Expected behavior
        }
    }

    @Test
    public void testWalkInOptimizedOrderChangingVisitor2() {
        final SparseFieldVector<Fraction> v = create(5);
        final FieldVectorChangingVisitor<Fraction> visitor;
        visitor = new FieldVectorChangingVisitor<Fraction>() {

            public Fraction visit(int index, Fraction value) {
                return Fraction.ZERO;
            }

            public void start(int dimension, int start, int end) {
                // Do nothing
            }

            public Fraction end() {
                return Fraction.ZERO;
            }
        };
        try {
            v.walkInOptimizedOrder(visitor, -1, 4);
            Assert.fail();
        } catch (OutOfRangeException e) {
            // Expected behavior
        }
        try {
            v.walkInOptimizedOrder(visitor, 5, 4);
            Assert.fail();
        } catch (OutOfRangeException e) {
            // Expected behavior
        }
        try {
            v.walkInOptimizedOrder(visitor, 0, -1);
            Assert.fail();
        } catch (OutOfRangeException e) {
            // Expected behavior
        }
        try {
            v.walkInOptimizedOrder(visitor, 0, 5);
            Assert.fail();
        } catch (OutOfRangeException e) {
            // Expected behavior
        }
        try {
            v.walkInOptimizedOrder(visitor, 4, 0);
            Assert.fail();
        } catch (NumberIsTooSmallException e) {
            // Expected behavior
        }
    }

[2962, 'src/test/java', 'org.apache.commons.math3.geometry.spherical.oned', 'ArcsSetTest', 'testArc', 39, 53]

[2962, 'src/test/java', 'org.apache.commons.math3.geometry.spherical.oned', 'ArcsSetTest', 'testWrapAround2PiArc', 55, 69]

    @Test
    public void testArc() {
        ArcsSet set = new ArcsSet(2.3, 5.7, 1.0e-10);
        Assert.assertEquals(3.4, set.getSize(), 1.0e-10);
        Assert.assertEquals(1.0e-10, set.getTolerance(), 1.0e-20);
        Assert.assertEquals(Region.Location.BOUNDARY, set.checkPoint(new S1Point(2.3)));
        Assert.assertEquals(Region.Location.BOUNDARY, set.checkPoint(new S1Point(5.7)));
        Assert.assertEquals(Region.Location.OUTSIDE,  set.checkPoint(new S1Point(1.2)));
        Assert.assertEquals(Region.Location.OUTSIDE,  set.checkPoint(new S1Point(8.5)));
        Assert.assertEquals(Region.Location.INSIDE,   set.checkPoint(new S1Point(8.7)));
        Assert.assertEquals(Region.Location.INSIDE,   set.checkPoint(new S1Point(3.0)));
        Assert.assertEquals(1, set.asList().size());
        Assert.assertEquals(2.3, set.asList().get(0).getInf(), 1.0e-10);
        Assert.assertEquals(5.7, set.asList().get(0).getSup(), 1.0e-10);
    }

    @Test
    public void testWrapAround2PiArc() {
        ArcsSet set = new ArcsSet(5.7 - MathUtils.TWO_PI, 2.3, 1.0e-10);
        Assert.assertEquals(MathUtils.TWO_PI - 3.4, set.getSize(), 1.0e-10);
        Assert.assertEquals(1.0e-10, set.getTolerance(), 1.0e-20);
        Assert.assertEquals(Region.Location.BOUNDARY, set.checkPoint(new S1Point(2.3)));
        Assert.assertEquals(Region.Location.BOUNDARY, set.checkPoint(new S1Point(5.7)));
        Assert.assertEquals(Region.Location.INSIDE,   set.checkPoint(new S1Point(1.2)));
        Assert.assertEquals(Region.Location.INSIDE,   set.checkPoint(new S1Point(8.5)));
        Assert.assertEquals(Region.Location.OUTSIDE,  set.checkPoint(new S1Point(8.7)));
        Assert.assertEquals(Region.Location.OUTSIDE,  set.checkPoint(new S1Point(3.0)));
        Assert.assertEquals(1, set.asList().size());
        Assert.assertEquals(5.7, set.asList().get(0).getInf(), 1.0e-10);
        Assert.assertEquals(2.3 + MathUtils.TWO_PI, set.asList().get(0).getSup(), 1.0e-10);
    }

[2999, 'src/test/java', 'org.apache.commons.math3.linear', 'RealVectorAbstractTest', 'testCosineLeftNullVector', 1315, 1320]

[2999, 'src/test/java', 'org.apache.commons.math3.linear', 'RealVectorAbstractTest', 'testCosineRightNullVector', 1322, 1327]

    @Test(expected=MathArithmeticException.class)
    public void testCosineLeftNullVector() {
        final RealVector v = create(new double[] {0, 0, 0});
        final RealVector w = create(new double[] {1, 0, 0});
        v.cosine(w);
    }

    @Test(expected=MathArithmeticException.class)
    public void testCosineRightNullVector() {
        final RealVector v = create(new double[] {0, 0, 0});
        final RealVector w = create(new double[] {1, 0, 0});
        w.cosine(v);
    }

[3019, 'src/test/java', 'org.apache.commons.math3.ode.nonstiff', 'ClassicalRungeKuttaIntegratorTest', 'testDecreasingSteps', 133, 182]

[3019, 'src/test/java', 'org.apache.commons.math3.ode.nonstiff', 'LutherIntegratorTest', 'testDecreasingSteps', 133, 182]

  @Test
  public void testDecreasingSteps()
      throws DimensionMismatchException, NumberIsTooSmallException,
             MaxCountExceededException, NoBracketingException {

    for (TestProblemAbstract pb : new TestProblemAbstract[] {
        new TestProblem1(), new TestProblem2(), new TestProblem3(),
        new TestProblem4(), new TestProblem5(), new TestProblem6()
    }) {

      double previousValueError = Double.NaN;
      double previousTimeError = Double.NaN;
      for (int i = 4; i < 10; ++i) {

        double step = (pb.getFinalTime() - pb.getInitialTime()) * FastMath.pow(2.0, -i);

        FirstOrderIntegrator integ = new ClassicalRungeKuttaIntegrator(step);
        TestProblemHandler handler = new TestProblemHandler(pb, integ);
        integ.addStepHandler(handler);
        EventHandler[] functions = pb.getEventsHandlers();
        for (int l = 0; l < functions.length; ++l) {
          integ.addEventHandler(functions[l],
                                     Double.POSITIVE_INFINITY, 1.0e-6 * step, 1000);
        }
        Assert.assertEquals(functions.length, integ.getEventHandlers().size());
        double stopTime = integ.integrate(pb, pb.getInitialTime(), pb.getInitialState(),
                                          pb.getFinalTime(), new double[pb.getDimension()]);
        if (functions.length == 0) {
            Assert.assertEquals(pb.getFinalTime(), stopTime, 1.0e-10);
        }

        double error = handler.getMaximalValueError();
        if (i > 4) {
          Assert.assertTrue(error < 1.01 * FastMath.abs(previousValueError));
        }
        previousValueError = error;

        double timeError = handler.getMaximalTimeError();
        if (i > 4) {
          Assert.assertTrue(timeError <= FastMath.abs(previousTimeError));
        }
        previousTimeError = timeError;

        integ.clearEventHandlers();
        Assert.assertEquals(0, integ.getEventHandlers().size());
      }

    }

  }

    @Test
    public void testDecreasingSteps()
            throws DimensionMismatchException, NumberIsTooSmallException,
            MaxCountExceededException, NoBracketingException {

        for (TestProblemAbstract pb : new TestProblemAbstract[] {
            new TestProblem1(), new TestProblem2(), new TestProblem3(),
            new TestProblem4(), new TestProblem5(), new TestProblem6()
        }) {

            double previousValueError = Double.NaN;
            double previousTimeError = Double.NaN;
            for (int i = 4; i < 10; ++i) {

                double step = (pb.getFinalTime() - pb.getInitialTime()) * FastMath.pow(2.0, -i);

                FirstOrderIntegrator integ = new LutherIntegrator(step);
                TestProblemHandler handler = new TestProblemHandler(pb, integ);
                integ.addStepHandler(handler);
                EventHandler[] functions = pb.getEventsHandlers();
                for (int l = 0; l < functions.length; ++l) {
                    integ.addEventHandler(functions[l],
                                          Double.POSITIVE_INFINITY, 1.0e-6 * step, 1000);
                }
                Assert.assertEquals(functions.length, integ.getEventHandlers().size());
                double stopTime = integ.integrate(pb, pb.getInitialTime(), pb.getInitialState(),
                                                  pb.getFinalTime(), new double[pb.getDimension()]);
                if (functions.length == 0) {
                    Assert.assertEquals(pb.getFinalTime(), stopTime, 1.0e-10);
                }

                double error = handler.getMaximalValueError();
                if (i > 4) {
                    Assert.assertTrue(error < 1.01 * FastMath.abs(previousValueError));
                }
                previousValueError = error;

                double timeError = handler.getMaximalTimeError();
                if (i > 4) {
                    Assert.assertTrue(timeError <= FastMath.abs(previousTimeError));
                }
                previousTimeError = timeError;

                integ.clearEventHandlers();
                Assert.assertEquals(0, integ.getEventHandlers().size());
            }

        }

    }

[3023, 'src/test/java', 'org.apache.commons.math3.analysis.polynomials', 'PolynomialFunctionTest', 'testfirstDerivativeComparison', 134, 154]

[3023, 'src/test/java', 'org.apache.commons.math3.analysis.polynomials', 'PolynomialFunctionTest', 'testMath341', 233, 253]

    @Test
    public void testfirstDerivativeComparison() {
        double[] f_coeff = { 3, 6, -2, 1 };
        double[] g_coeff = { 6, -4, 3 };
        double[] h_coeff = { -4, 6 };

        PolynomialFunction f = new PolynomialFunction(f_coeff);
        PolynomialFunction g = new PolynomialFunction(g_coeff);
        PolynomialFunction h = new PolynomialFunction(h_coeff);

        // compare f' = g
        Assert.assertEquals(f.derivative().value(0), g.value(0), tolerance);
        Assert.assertEquals(f.derivative().value(1), g.value(1), tolerance);
        Assert.assertEquals(f.derivative().value(100), g.value(100), tolerance);
        Assert.assertEquals(f.derivative().value(4.1), g.value(4.1), tolerance);
        Assert.assertEquals(f.derivative().value(-3.25), g.value(-3.25), tolerance);

        // compare g' = h
        Assert.assertEquals(g.derivative().value(FastMath.PI), h.value(FastMath.PI), tolerance);
        Assert.assertEquals(g.derivative().value(FastMath.E),  h.value(FastMath.E),  tolerance);
    }

    @Test
    public void testMath341() {
        double[] f_coeff = { 3, 6, -2, 1 };
        double[] g_coeff = { 6, -4, 3 };
        double[] h_coeff = { -4, 6 };

        PolynomialFunction f = new PolynomialFunction(f_coeff);
        PolynomialFunction g = new PolynomialFunction(g_coeff);
        PolynomialFunction h = new PolynomialFunction(h_coeff);

        // compare f' = g
        Assert.assertEquals(f.derivative().value(0), g.value(0), tolerance);
        Assert.assertEquals(f.derivative().value(1), g.value(1), tolerance);
        Assert.assertEquals(f.derivative().value(100), g.value(100), tolerance);
        Assert.assertEquals(f.derivative().value(4.1), g.value(4.1), tolerance);
        Assert.assertEquals(f.derivative().value(-3.25), g.value(-3.25), tolerance);

        // compare g' = h
        Assert.assertEquals(g.derivative().value(FastMath.PI), h.value(FastMath.PI), tolerance);
        Assert.assertEquals(g.derivative().value(FastMath.E),  h.value(FastMath.E),  tolerance);
    }

[3037, 'src/test/java', 'org.apache.commons.math3.fraction', 'BigFractionFormatTest', 'testParseProperInvalidMinus', 240, 256]

[3037, 'src/test/java', 'org.apache.commons.math3.fraction', 'FractionFormatTest', 'testParseProperInvalidMinus', 278, 294]

    @Test
    public void testParseProperInvalidMinus() {
        String source = "2 -2 / 3";
        try {
            properFormat.parse(source);
            Assert.fail("invalid minus in improper fraction.");
        } catch (MathParseException ex) {
            // expected
        }
        source = "2 2 / -3";
        try {
            properFormat.parse(source);
            Assert.fail("invalid minus in improper fraction.");
        } catch (MathParseException ex) {
            // expected
        }
    }

    @Test
    public void testParseProperInvalidMinus() {
        String source = "2 -2 / 3";
        try {
            properFormat.parse(source);
            Assert.fail("invalid minus in improper fraction.");
        } catch (MathParseException ex) {
            // expected
        }
        source = "2 2 / -3";
        try {
            properFormat.parse(source);
            Assert.fail("invalid minus in improper fraction.");
        } catch (MathParseException ex) {
            // expected
        }
    }

[3038, 'src/test/java', 'org.apache.commons.math3.complex', 'ComplexTest', 'testEqualsRealDifference', 604, 609]

[3038, 'src/test/java', 'org.apache.commons.math3.complex', 'ComplexTest', 'testEqualsImaginaryDifference', 611, 616]

    @Test
    public void testEqualsRealDifference() {
        Complex x = new Complex(0.0, 0.0);
        Complex y = new Complex(0.0 + Double.MIN_VALUE, 0.0);
        Assert.assertFalse(x.equals(y));
    }

    @Test
    public void testEqualsImaginaryDifference() {
        Complex x = new Complex(0.0, 0.0);
        Complex y = new Complex(0.0, 0.0 + Double.MIN_VALUE);
        Assert.assertFalse(x.equals(y));
    }

[3045, 'src/test/java', 'org.apache.commons.math3.transform', 'FastFourierTransformerTest', 'testTransformComplex', 335, 352]

[3045, 'src/test/java', 'org.apache.commons.math3.transform', 'FastFourierTransformerTest', 'testStandardTransformReal', 354, 371]

    @Test
    public void testTransformComplex() {
        final DftNormalization[] norm;
        norm = DftNormalization.values();
        final TransformType[] type;
        type = TransformType.values();
        for (int i = 0; i < norm.length; i++) {
            for (int j = 0; j < type.length; j++) {
                doTestTransformComplex(2, 1.0E-15, norm[i], type[j]);
                doTestTransformComplex(4, 1.0E-14, norm[i], type[j]);
                doTestTransformComplex(8, 1.0E-14, norm[i], type[j]);
                doTestTransformComplex(16, 1.0E-13, norm[i], type[j]);
                doTestTransformComplex(32, 1.0E-13, norm[i], type[j]);
                doTestTransformComplex(64, 1.0E-12, norm[i], type[j]);
                doTestTransformComplex(128, 1.0E-12, norm[i], type[j]);
            }
        }
    }

    @Test
    public void testStandardTransformReal() {
        final DftNormalization[] norm;
        norm = DftNormalization.values();
        final TransformType[] type;
        type = TransformType.values();
        for (int i = 0; i < norm.length; i++) {
            for (int j = 0; j < type.length; j++) {
                doTestTransformReal(2, 1.0E-15, norm[i], type[j]);
                doTestTransformReal(4, 1.0E-14, norm[i], type[j]);
                doTestTransformReal(8, 1.0E-14, norm[i], type[j]);
                doTestTransformReal(16, 1.0E-13, norm[i], type[j]);
                doTestTransformReal(32, 1.0E-13, norm[i], type[j]);
                doTestTransformReal(64, 1.0E-13, norm[i], type[j]);
                doTestTransformReal(128, 1.0E-11, norm[i], type[j]);
            }
        }
    }

[3048, 'src/test/java', 'org.apache.commons.math3.util', 'PrecisionTest', 'testEqualsWithAllowedDelta', 77, 88]

[3048, 'src/test/java', 'org.apache.commons.math3.util', 'PrecisionTest', 'testEqualsIncludingNaNWithAllowedDelta', 109, 120]

    @Test
    public void testEqualsWithAllowedDelta() {
        Assert.assertTrue(Precision.equals(153.0000, 153.0000, .0625));
        Assert.assertTrue(Precision.equals(153.0000, 153.0625, .0625));
        Assert.assertTrue(Precision.equals(152.9375, 153.0000, .0625));
        Assert.assertFalse(Precision.equals(153.0000, 153.0625, .0624));
        Assert.assertFalse(Precision.equals(152.9374, 153.0000, .0625));
        Assert.assertFalse(Precision.equals(Double.NaN, Double.NaN, 1.0));
        Assert.assertTrue(Precision.equals(Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY, 1.0));
        Assert.assertTrue(Precision.equals(Double.NEGATIVE_INFINITY, Double.NEGATIVE_INFINITY, 1.0));
        Assert.assertFalse(Precision.equals(Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY, 1.0));
    }

    @Test
    public void testEqualsIncludingNaNWithAllowedDelta() {
        Assert.assertTrue(Precision.equalsIncludingNaN(153.0000, 153.0000, .0625));
        Assert.assertTrue(Precision.equalsIncludingNaN(153.0000, 153.0625, .0625));
        Assert.assertTrue(Precision.equalsIncludingNaN(152.9375, 153.0000, .0625));
        Assert.assertTrue(Precision.equalsIncludingNaN(Double.NaN, Double.NaN, 1.0));
        Assert.assertTrue(Precision.equalsIncludingNaN(Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY, 1.0));
        Assert.assertTrue(Precision.equalsIncludingNaN(Double.NEGATIVE_INFINITY, Double.NEGATIVE_INFINITY, 1.0));
        Assert.assertFalse(Precision.equalsIncludingNaN(Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY, 1.0));
        Assert.assertFalse(Precision.equalsIncludingNaN(153.0000, 153.0625, .0624));
        Assert.assertFalse(Precision.equalsIncludingNaN(152.9374, 153.0000, .0625));
    }

[3056, 'src/test/java', 'org.apache.commons.math3.analysis.integration', 'IterativeLegendreGaussIntegratorTest', 'testExactIntegration', 78, 104]

[3056, 'src/test/java', 'org.apache.commons.math3.analysis.integration', 'LegendreGaussIntegratorTest', 'testExactIntegration', 77, 103]

    @Test
    public void testExactIntegration() {
        Random random = new Random(86343623467878363l);
        for (int n = 2; n < 6; ++n) {
            IterativeLegendreGaussIntegrator integrator =
                new IterativeLegendreGaussIntegrator(n,
                                                     BaseAbstractUnivariateIntegrator.DEFAULT_RELATIVE_ACCURACY,
                                                     BaseAbstractUnivariateIntegrator.DEFAULT_ABSOLUTE_ACCURACY,
                                                     BaseAbstractUnivariateIntegrator.DEFAULT_MIN_ITERATIONS_COUNT,
                                                     64);

            // an n points Gauss-Legendre integrator integrates 2n-1 degree polynoms exactly
            for (int degree = 0; degree <= 2 * n - 1; ++degree) {
                for (int i = 0; i < 10; ++i) {
                    double[] coeff = new double[degree + 1];
                    for (int k = 0; k < coeff.length; ++k) {
                        coeff[k] = 2 * random.nextDouble() - 1;
                    }
                    PolynomialFunction p = new PolynomialFunction(coeff);
                    double result    = integrator.integrate(10000, p, -5.0, 15.0);
                    double reference = exactIntegration(p, -5.0, 15.0);
                    Assert.assertEquals(n + " " + degree + " " + i, reference, result, 1.0e-12 * (1.0 + FastMath.abs(reference)));
                }
            }

        }
    }

    @Test
    public void testExactIntegration() {
        Random random = new Random(86343623467878363l);
        for (int n = 2; n < 6; ++n) {
            LegendreGaussIntegrator integrator =
                new LegendreGaussIntegrator(n,
                                            BaseAbstractUnivariateIntegrator.DEFAULT_RELATIVE_ACCURACY,
                                            BaseAbstractUnivariateIntegrator.DEFAULT_ABSOLUTE_ACCURACY,
                                            BaseAbstractUnivariateIntegrator.DEFAULT_MIN_ITERATIONS_COUNT,
                                            64);

            // an n points Gauss-Legendre integrator integrates 2n-1 degree polynoms exactly
            for (int degree = 0; degree <= 2 * n - 1; ++degree) {
                for (int i = 0; i < 10; ++i) {
                    double[] coeff = new double[degree + 1];
                    for (int k = 0; k < coeff.length; ++k) {
                        coeff[k] = 2 * random.nextDouble() - 1;
                    }
                    PolynomialFunction p = new PolynomialFunction(coeff);
                    double result    = integrator.integrate(10000, p, -5.0, 15.0);
                    double reference = exactIntegration(p, -5.0, 15.0);
                    Assert.assertEquals(n + " " + degree + " " + i, reference, result, 1.0e-12 * (1.0 + FastMath.abs(reference)));
                }
            }

        }
    }

[3086, 'src/test/java', 'org.apache.commons.math3.linear', 'BlockRealMatrixTest', 'testGetSubMatrix', 497, 515]

[3086, 'src/test/java', 'org.apache.commons.math3.linear', 'BlockRealMatrixTest', 'testCopySubMatrix', 585, 604]

    @Test
    public void testGetSubMatrix() {
        RealMatrix m = new BlockRealMatrix(subTestData);
        checkGetSubMatrix(m, subRows23Cols00,  2 , 3 , 0, 0);
        checkGetSubMatrix(m, subRows00Cols33,  0 , 0 , 3, 3);
        checkGetSubMatrix(m, subRows01Cols23,  0 , 1 , 2, 3);
        checkGetSubMatrix(m, subRows02Cols13,  new int[] { 0, 2 }, new int[] { 1, 3 });
        checkGetSubMatrix(m, subRows03Cols12,  new int[] { 0, 3 }, new int[] { 1, 2 });
        checkGetSubMatrix(m, subRows03Cols123, new int[] { 0, 3 }, new int[] { 1, 2, 3 });
        checkGetSubMatrix(m, subRows20Cols123, new int[] { 2, 0 }, new int[] { 1, 2, 3 });
        checkGetSubMatrix(m, subRows31Cols31,  new int[] { 3, 1 }, new int[] { 3, 1 });
        checkGetSubMatrix(m, subRows31Cols31,  new int[] { 3, 1 }, new int[] { 3, 1 });
        checkGetSubMatrix(m, null,  1, 0, 2, 4);
        checkGetSubMatrix(m, null, -1, 1, 2, 2);
        checkGetSubMatrix(m, null,  1, 0, 2, 2);
        checkGetSubMatrix(m, null,  1, 0, 2, 4);
        checkGetSubMatrix(m, null, new int[] {},    new int[] { 0 });
        checkGetSubMatrix(m, null, new int[] { 0 }, new int[] { 4 });
    }

    @Test
    public void testCopySubMatrix() {
        RealMatrix m = new BlockRealMatrix(subTestData);
        checkCopy(m, subRows23Cols00,  2 , 3 , 0, 0);
        checkCopy(m, subRows00Cols33,  0 , 0 , 3, 3);
        checkCopy(m, subRows01Cols23,  0 , 1 , 2, 3);
        checkCopy(m, subRows02Cols13,  new int[] { 0, 2 }, new int[] { 1, 3 });
        checkCopy(m, subRows03Cols12,  new int[] { 0, 3 }, new int[] { 1, 2 });
        checkCopy(m, subRows03Cols123, new int[] { 0, 3 }, new int[] { 1, 2, 3 });
        checkCopy(m, subRows20Cols123, new int[] { 2, 0 }, new int[] { 1, 2, 3 });
        checkCopy(m, subRows31Cols31,  new int[] { 3, 1 }, new int[] { 3, 1 });
        checkCopy(m, subRows31Cols31,  new int[] { 3, 1 }, new int[] { 3, 1 });

        checkCopy(m, null,  1, 0, 2, 4);
        checkCopy(m, null, -1, 1, 2, 2);
        checkCopy(m, null,  1, 0, 2, 2);
        checkCopy(m, null,  1, 0, 2, 4);
        checkCopy(m, null, new int[] {},    new int[] { 0 });
        checkCopy(m, null, new int[] { 0 }, new int[] { 4 });
    }

[3112, 'src/test/java', 'org.apache.commons.math3.distribution', 'WeibullDistributionTest', 'testAlpha', 83, 93]

[3112, 'src/test/java', 'org.apache.commons.math3.distribution', 'WeibullDistributionTest', 'testBeta', 95, 105]

    @Test
    public void testAlpha() {
        WeibullDistribution dist = new WeibullDistribution(1, 2);
        Assert.assertEquals(1, dist.getShape(), 0);
        try {
            new WeibullDistribution(0, 2);
            Assert.fail("NotStrictlyPositiveException expected");
        } catch (NotStrictlyPositiveException e) {
            // Expected.
        }
    }

    @Test
    public void testBeta() {
        WeibullDistribution dist = new WeibullDistribution(1, 2);
        Assert.assertEquals(2, dist.getScale(), 0);
        try {
            new WeibullDistribution(1, 0);
            Assert.fail("NotStrictlyPositiveException expected");
        } catch (NotStrictlyPositiveException e) {
            // Expected.
        }
    }

[3117, 'src/test/java', 'org.apache.commons.math3.linear', 'HessenbergTransformerTest', 'testRandomData', 86, 105]

[3117, 'src/test/java', 'org.apache.commons.math3.linear', 'SchurTransformerTest', 'testRandomData', 90, 109]

    @Test
    public void testRandomData() {
        for (int run = 0; run < 100; run++) {
            Random r = new Random(System.currentTimeMillis());

            // matrix size
            int size = r.nextInt(20) + 4;

            double[][] data = new double[size][size];
            for (int i = 0; i < size; i++) {
                for (int j = 0; j < size; j++) {
                    data[i][j] = r.nextInt(100);
                }
            }

            RealMatrix m = MatrixUtils.createRealMatrix(data);
            RealMatrix h = checkAEqualPHPt(m);
            checkHessenbergForm(h);
        }
    }

    @Test
    public void testRandomData() {
        for (int run = 0; run < 100; run++) {
            Random r = new Random(System.currentTimeMillis());

            // matrix size
            int size = r.nextInt(20) + 4;

            double[][] data = new double[size][size];
            for (int i = 0; i < size; i++) {
                for (int j = 0; j < size; j++) {
                    data[i][j] = r.nextInt(100);
                }
            }

            RealMatrix m = MatrixUtils.createRealMatrix(data);
            RealMatrix s = checkAEqualPTPt(m);
            checkSchurForm(s);
        }
    }

[3129, 'src/test/java', 'org.apache.commons.math3.stat.inference', 'ChiSquareTestTest', 'testChiSquareDataSetsComparisonBadCounts', 229, 259]

[3129, 'src/test/java', 'org.apache.commons.math3.stat.inference', 'GTestTest', 'testGTestSetsComparisonBadCounts', 153, 182]

    @Test
    public void testChiSquareDataSetsComparisonBadCounts()
        {
        long[] observed1 = {10, -1, 12, 10, 15};
        long[] observed2 = {15, 10, 10, 15, 5};
        try {
            testStatistic.chiSquareTestDataSetsComparison(
                    observed1, observed2);
            Assert.fail("Expecting NotPositiveException - negative count");
        } catch (NotPositiveException ex) {
            // expected
        }
        long[] observed3 = {10, 0, 12, 10, 15};
        long[] observed4 = {15, 0, 10, 15, 5};
        try {
            testStatistic.chiSquareTestDataSetsComparison(
                    observed3, observed4);
            Assert.fail("Expecting ZeroException - double 0's");
        } catch (ZeroException ex) {
            // expected
        }
        long[] observed5 = {10, 10, 12, 10, 15};
        long[] observed6 = {0, 0, 0, 0, 0};
        try {
            testStatistic.chiSquareTestDataSetsComparison(
                    observed5, observed6);
            Assert.fail("Expecting ZeroException - vanishing counts");
        } catch (ZeroException ex) {
            // expected
        }
    }

    @Test
    public void testGTestSetsComparisonBadCounts() {
        long[] observed1 = {10, -1, 12, 10, 15};
        long[] observed2 = {15, 10, 10, 15, 5};
        try {
            testStatistic.gTestDataSetsComparison(
                    observed1, observed2);
            Assert.fail("Expecting NotPositiveException - negative count");
        } catch (NotPositiveException ex) {
            // expected
        }
        long[] observed3 = {10, 0, 12, 10, 15};
        long[] observed4 = {15, 0, 10, 15, 5};
        try {
            testStatistic.gTestDataSetsComparison(
                    observed3, observed4);
            Assert.fail("Expecting ZeroException - double 0's");
        } catch (ZeroException ex) {
            // expected
        }
        long[] observed5 = {10, 10, 12, 10, 15};
        long[] observed6 = {0, 0, 0, 0, 0};
        try {
            testStatistic.gTestDataSetsComparison(
                    observed5, observed6);
            Assert.fail("Expecting ZeroException - vanishing counts");
        } catch (ZeroException ex) {
            // expected
        }
    }

[3161, 'src/test/java', 'org.apache.commons.math3.linear', 'RealMatrixFormatAbstractTest', 'testParseNegativeComponent', 238, 248]

[3161, 'src/test/java', 'org.apache.commons.math3.linear', 'RealMatrixFormatAbstractTest', 'testParseZeroComponent', 262, 272]

    @Test
    public void testParseNegativeComponent() {
        String source =
            "{{-1" + getDecimalCharacter() +
            "2323,1" + getDecimalCharacter() +
            "4343,1" + getDecimalCharacter() +
            "6333}}";
        RealMatrix expected = MatrixUtils.createRealMatrix(new double[][] {{-1.2323, 1.4343, 1.6333}});
        RealMatrix actual = realMatrixFormat.parse(source);
        Assert.assertEquals(expected, actual);
    }

    @Test
    public void testParseZeroComponent() {
        String source =
            "{{0" + getDecimalCharacter() +
            "0,-1" + getDecimalCharacter() +
            "4343,1" + getDecimalCharacter() +
            "6333}}";
        RealMatrix expected = MatrixUtils.createRealMatrix(new double[][] {{0.0, -1.4343, 1.6333}});
        RealMatrix actual = realMatrixFormat.parse(source);
        Assert.assertEquals(expected, actual);
    }

[3163, 'src/test/java', 'org.apache.commons.math3.stat.inference', 'TTestTest', 'testSmallSamples', 272, 282]

[3163, 'src/test/java', 'org.apache.commons.math3.stat.inference', 'TestUtilsTest', 'testSmallSamples', 425, 435]

    @Test
    public void testSmallSamples() {
        double[] sample1 = {1d, 3d};
        double[] sample2 = {4d, 5d};

        // Target values computed using R, version 1.8.1 (linux version)
        Assert.assertEquals(-2.2360679775, testStatistic.t(sample1, sample2),
                1E-10);
        Assert.assertEquals(0.198727388935, testStatistic.tTest(sample1, sample2),
                1E-10);
    }

    @Test
    public void testSmallSamples() {
        double[] sample1 = {1d, 3d};
        double[] sample2 = {4d, 5d};

        // Target values computed using R, version 1.8.1 (linux version)
        Assert.assertEquals(-2.2360679775, TestUtils.t(sample1, sample2),
                1E-10);
        Assert.assertEquals(0.198727388935, TestUtils.tTest(sample1, sample2),
                1E-10);
    }

[3174, 'src/test/java', 'org.apache.commons.math3.distribution', 'BinomialDistributionTest', 'testDegenerate1', 112, 127]

[3174, 'src/test/java', 'org.apache.commons.math3.distribution', 'BinomialDistributionTest', 'testDegenerate2', 130, 145]

    @Test
    public void testDegenerate1() {
        BinomialDistribution dist = new BinomialDistribution(5, 1.0d);
        setDistribution(dist);
        setCumulativeTestPoints(new int[] { -1, 0, 1, 2, 5, 10 });
        setCumulativeTestValues(new double[] { 0d, 0d, 0d, 0d, 1d, 1d });
        setDensityTestPoints(new int[] { -1, 0, 1, 2, 5, 10 });
        setDensityTestValues(new double[] { 0d, 0d, 0d, 0d, 1d, 0d });
        setInverseCumulativeTestPoints(new double[] { 0.1d, 0.5d });
        setInverseCumulativeTestValues(new int[] { 5, 5 });
        verifyDensities();
        verifyCumulativeProbabilities();
        verifyInverseCumulativeProbabilities();
        Assert.assertEquals(dist.getSupportLowerBound(), 5);
        Assert.assertEquals(dist.getSupportUpperBound(), 5);
    }

    @Test
    public void testDegenerate2() {
        BinomialDistribution dist = new BinomialDistribution(0, 0.01d);
        setDistribution(dist);
        setCumulativeTestPoints(new int[] { -1, 0, 1, 2, 5, 10 });
        setCumulativeTestValues(new double[] { 0d, 1d, 1d, 1d, 1d, 1d });
        setDensityTestPoints(new int[] { -1, 0, 1, 2, 5, 10 });
        setDensityTestValues(new double[] { 0d, 1d, 0d, 0d, 0d, 0d });
        setInverseCumulativeTestPoints(new double[] { 0.1d, 0.5d });
        setInverseCumulativeTestValues(new int[] { 0, 0 });
        verifyDensities();
        verifyCumulativeProbabilities();
        verifyInverseCumulativeProbabilities();
        Assert.assertEquals(dist.getSupportLowerBound(), 0);
        Assert.assertEquals(dist.getSupportUpperBound(), 0);
    }

[3244, 'src/test/java', 'org.apache.commons.math3.linear', 'Array2DRowRealMatrixTest', 'testOperate', 302, 316]

[3244, 'src/test/java', 'org.apache.commons.math3.linear', 'Array2DRowRealMatrixTest', 'testPremultiplyVector', 344, 358]

    @Test
    public void testOperate() {
        RealMatrix m = new Array2DRowRealMatrix(id);
        TestUtils.assertEquals("identity operate", testVector,
                    m.operate(testVector), entryTolerance);
        TestUtils.assertEquals("identity operate", testVector,
                    m.operate(new ArrayRealVector(testVector)).toArray(), entryTolerance);
        m = new Array2DRowRealMatrix(bigSingular);
        try {
            m.operate(testVector);
            Assert.fail("Expecting illegalArgumentException");
        } catch (MathIllegalArgumentException ex) {
            // ignored
        }
    }

    @Test
    public void testPremultiplyVector() {
        RealMatrix m = new Array2DRowRealMatrix(testData);
        TestUtils.assertEquals("premultiply", m.preMultiply(testVector),
                    preMultTest, normTolerance);
        TestUtils.assertEquals("premultiply", m.preMultiply(new ArrayRealVector(testVector).toArray()),
                    preMultTest, normTolerance);
        m = new Array2DRowRealMatrix(bigSingular);
        try {
            m.preMultiply(testVector);
            Assert.fail("expecting MathIllegalArgumentException");
        } catch (MathIllegalArgumentException ex) {
            // ignored
        }
    }

[3248, 'src/test/java', 'org.apache.commons.math3.fraction', 'BigFractionTest', 'testFloatValue', 241, 248]

[3248, 'src/test/java', 'org.apache.commons.math3.fraction', 'FractionTest', 'testFloatValue', 197, 204]

    @Test
    public void testFloatValue() {
        BigFraction first = new BigFraction(1, 2);
        BigFraction second = new BigFraction(1, 3);

        Assert.assertEquals(0.5f, first.floatValue(), 0.0f);
        Assert.assertEquals((float) (1.0 / 3.0), second.floatValue(), 0.0f);
    }

    @Test
    public void testFloatValue() {
        Fraction first = new Fraction(1, 2);
        Fraction second = new Fraction(1, 3);

        Assert.assertEquals(0.5f, first.floatValue(), 0.0f);
        Assert.assertEquals((float)(1.0 / 3.0), second.floatValue(), 0.0f);
    }

[3279, 'src/test/java', 'org.apache.commons.math3.fraction', 'BigFractionFormatTest', 'testParseProperNegative', 222, 238]

[3279, 'src/test/java', 'org.apache.commons.math3.fraction', 'FractionFormatTest', 'testParseProperNegative', 260, 276]

    @Test
    public void testParseProperNegative() {
        String source = "-1 2 / 3";
        {
            BigFraction c = properFormat.parse(source);
            Assert.assertNotNull(c);
            Assert.assertEquals(-5, c.getNumeratorAsInt());
            Assert.assertEquals(3, c.getDenominatorAsInt());
        }

        try {
            improperFormat.parse(source);
            Assert.fail("invalid improper fraction.");
        } catch (MathParseException ex) {
            // success
        }
    }

    @Test
    public void testParseProperNegative() {
        String source = "-1 2 / 3";
        {
            Fraction c = properFormat.parse(source);
            Assert.assertNotNull(c);
            Assert.assertEquals(-5, c.getNumerator());
            Assert.assertEquals(3, c.getDenominator());
        }

        try {
            improperFormat.parse(source);
            Assert.fail("invalid improper fraction.");
        } catch (MathParseException ex) {
            // success
        }
    }

[3281, 'src/test/java', 'org.apache.commons.math3.analysis.function', 'GaussianTest', 'testParametricValue', 131, 142]

[3281, 'src/test/java', 'org.apache.commons.math3.analysis.function', 'HarmonicOscillatorTest', 'testParametricValue', 109, 120]

    @Test
    public void testParametricValue() {
        final double norm = 2;
        final double mean = 3;
        final double sigma = 4;
        final Gaussian f = new Gaussian(norm, mean, sigma);

        final Gaussian.Parametric g = new Gaussian.Parametric();
        Assert.assertEquals(f.value(-1), g.value(-1, new double[] {norm, mean, sigma}), 0);
        Assert.assertEquals(f.value(0), g.value(0, new double[] {norm, mean, sigma}), 0);
        Assert.assertEquals(f.value(2), g.value(2, new double[] {norm, mean, sigma}), 0);
    }

    @Test
    public void testParametricValue() {
        final double amplitude = 2;
        final double omega = 3;
        final double phase = 4;
        final HarmonicOscillator f = new HarmonicOscillator(amplitude, omega, phase);

        final HarmonicOscillator.Parametric g = new HarmonicOscillator.Parametric();
        Assert.assertEquals(f.value(-1), g.value(-1, new double[] {amplitude, omega, phase}), 0);
        Assert.assertEquals(f.value(0), g.value(0, new double[] {amplitude, omega, phase}), 0);
        Assert.assertEquals(f.value(2), g.value(2, new double[] {amplitude, omega, phase}), 0);
    }

[3282, 'src/test/java', 'org.apache.commons.math3.stat.descriptive.rank', 'MaxTest', 'testSpecialValues', 47, 62]

[3282, 'src/test/java', 'org.apache.commons.math3.stat.descriptive.rank', 'MinTest', 'testSpecialValues', 47, 62]

    @Test
    public void testSpecialValues() {
        double[] testArray = {0d, Double.NaN, Double.NEGATIVE_INFINITY,
                Double.POSITIVE_INFINITY};
        Max max = new Max();
        Assert.assertTrue(Double.isNaN(max.getResult()));
        max.increment(testArray[0]);
        Assert.assertEquals(0d, max.getResult(), 0);
        max.increment(testArray[1]);
        Assert.assertEquals(0d, max.getResult(), 0);
        max.increment(testArray[2]);
        Assert.assertEquals(0d, max.getResult(), 0);
        max.increment(testArray[3]);
        Assert.assertEquals(Double.POSITIVE_INFINITY, max.getResult(), 0);
        Assert.assertEquals(Double.POSITIVE_INFINITY, max.evaluate(testArray), 0);
    }

    @Test
    public void testSpecialValues() {
        double[] testArray = {0d, Double.NaN, Double.POSITIVE_INFINITY,
                Double.NEGATIVE_INFINITY};
        Min min = new Min();
        Assert.assertTrue(Double.isNaN(min.getResult()));
        min.increment(testArray[0]);
        Assert.assertEquals(0d, min.getResult(), 0);
        min.increment(testArray[1]);
        Assert.assertEquals(0d, min.getResult(), 0);
        min.increment(testArray[2]);
        Assert.assertEquals(0d, min.getResult(), 0);
        min.increment(testArray[3]);
        Assert.assertEquals(Double.NEGATIVE_INFINITY, min.getResult(), 0);
        Assert.assertEquals(Double.NEGATIVE_INFINITY, min.evaluate(testArray), 0);
    }

[3283, 'src/test/java', 'org.apache.commons.math3.stat.correlation', 'KendallsCorrelationTest', 'testSwiss', 66, 79]

[3283, 'src/test/java', 'org.apache.commons.math3.stat.correlation', 'SpearmansRankCorrelationTest', 'testSwiss', 60, 73]

    @Test
    public void testSwiss() {
        RealMatrix matrix = createRealMatrix(swissData, 47, 5);
        KendallsCorrelation corrInstance = new KendallsCorrelation(matrix);
        RealMatrix correlationMatrix = corrInstance.getCorrelationMatrix();
        double[] rData = new double[] {
                1, 0.1795465254708308, -0.4762437404200669, -0.3306111613580587, 0.2453703703703704,
                0.1795465254708308, 1, -0.4505221560842292, -0.4761645631778491, 0.2054604569820847,
                -0.4762437404200669, -0.4505221560842292, 1, 0.528943683925829, -0.3212755391722673,
                -0.3306111613580587, -0.4761645631778491, 0.528943683925829, 1, -0.08479652265379604,
                0.2453703703703704, 0.2054604569820847, -0.3212755391722673, -0.08479652265379604, 1
        };
        TestUtils.assertEquals("Kendall's correlation matrix", createRealMatrix(rData, 5, 5), correlationMatrix, 10E-15);
    }

    @Test
    public void testSwiss() {
        RealMatrix matrix = createRealMatrix(swissData, 47, 5);
        SpearmansCorrelation corrInstance = new SpearmansCorrelation(matrix);
        RealMatrix correlationMatrix = corrInstance.getCorrelationMatrix();
        double[] rData = new double[] {
                1, 0.2426642769364176, -0.660902996352354, -0.443257690360988, 0.4136455623012432,
                0.2426642769364176, 1, -0.598859938748963, -0.650463814145816, 0.2886878090882852,
               -0.660902996352354, -0.598859938748963, 1, 0.674603831406147, -0.4750575257171745,
               -0.443257690360988, -0.650463814145816, 0.674603831406147, 1, -0.1444163088302244,
                0.4136455623012432, 0.2886878090882852, -0.4750575257171745, -0.1444163088302244, 1
        };
        TestUtils.assertEquals("Spearman's correlation matrix", createRealMatrix(rData, 5, 5), correlationMatrix, 10E-15);
    }

[3294, 'src/test/java', 'org.apache.commons.math3.random', 'RandomDataGeneratorTest', 'testNextBeta', 1012, 1022]

[3294, 'src/test/java', 'org.apache.commons.math3.random', 'RandomDataGeneratorTest', 'testNextF', 1048, 1058]

    @Test
    public void testNextBeta() {
        double[] quartiles = TestUtils.getDistributionQuartiles(new BetaDistribution(2,5));
        long[] counts = new long[4];
        randomData.reSeed(1000);
        for (int i = 0; i < 1000; i++) {
            double value = randomData.nextBeta(2, 5);
            TestUtils.updateCounts(value, counts, quartiles);
        }
        TestUtils.assertChiSquareAccept(expected, counts, 0.001);
    }

    @Test
    public void testNextF() {
        double[] quartiles = TestUtils.getDistributionQuartiles(new FDistribution(12, 5));
        long[] counts = new long[4];
        randomData.reSeed(1000);
        for (int i = 0; i < 1000; i++) {
            double value = randomData.nextF(12, 5);
            TestUtils.updateCounts(value, counts, quartiles);
        }
        TestUtils.assertChiSquareAccept(expected, counts, 0.001);
    }

[3295, 'src/test/java', 'org.apache.commons.math3.random', 'RandomDataGeneratorTest', 'testNextCauchy', 1024, 1034]

[3295, 'src/test/java', 'org.apache.commons.math3.random', 'RandomDataGeneratorTest', 'testNextWeibull', 1098, 1108]

    @Test
    public void testNextCauchy() {
        double[] quartiles = TestUtils.getDistributionQuartiles(new CauchyDistribution(1.2, 2.1));
        long[] counts = new long[4];
        randomData.reSeed(1000);
        for (int i = 0; i < 1000; i++) {
            double value = randomData.nextCauchy(1.2, 2.1);
            TestUtils.updateCounts(value, counts, quartiles);
        }
        TestUtils.assertChiSquareAccept(expected, counts, 0.001);
    }

    @Test
    public void testNextWeibull() {
        double[] quartiles = TestUtils.getDistributionQuartiles(new WeibullDistribution(1.2, 2.1));
        long[] counts = new long[4];
        randomData.reSeed(1000);
        for (int i = 0; i < 1000; i++) {
            double value = randomData.nextWeibull(1.2, 2.1);
            TestUtils.updateCounts(value, counts, quartiles);
        }
        TestUtils.assertChiSquareAccept(expected, counts, 0.001);
    }

[3308, 'src/test/java', 'org.apache.commons.math3.linear', 'BlockFieldMatrixTest', 'testTranspose', 433, 442]

[3308, 'src/test/java', 'org.apache.commons.math3.linear', 'FieldMatrixImplTest', 'testTranspose', 296, 305]

    @Test
    public void testTranspose() {
        FieldMatrix<Fraction> m = new BlockFieldMatrix<Fraction>(testData);
        FieldMatrix<Fraction> mIT = new FieldLUDecomposition<Fraction>(m).getSolver().getInverse().transpose();
        FieldMatrix<Fraction> mTI = new FieldLUDecomposition<Fraction>(m.transpose()).getSolver().getInverse();
        TestUtils.assertEquals(mIT, mTI);
        m = new BlockFieldMatrix<Fraction>(testData2);
        FieldMatrix<Fraction> mt = new BlockFieldMatrix<Fraction>(testData2T);
        TestUtils.assertEquals(mt, m.transpose());
    }

    @Test
    public void testTranspose() {
        FieldMatrix<Fraction> m = new Array2DRowFieldMatrix<Fraction>(testData);
        FieldMatrix<Fraction> mIT = new FieldLUDecomposition<Fraction>(m).getSolver().getInverse().transpose();
        FieldMatrix<Fraction> mTI = new FieldLUDecomposition<Fraction>(m.transpose()).getSolver().getInverse();
        TestUtils.assertEquals(mIT, mTI);
        m = new Array2DRowFieldMatrix<Fraction>(testData2);
        FieldMatrix<Fraction> mt = new Array2DRowFieldMatrix<Fraction>(testData2T);
        TestUtils.assertEquals(mt, m.transpose());
    }

[3309, 'src/test/java', 'org.apache.commons.math3.random', 'HaltonSequenceGeneratorTest', 'testSkip', 122, 133]

[3309, 'src/test/java', 'org.apache.commons.math3.random', 'SobolSequenceGeneratorTest', 'testSkip', 94, 105]

    @Test
    public void testSkip() {
        double[] result = generator.skipTo(5);
        Assert.assertArrayEquals(referenceValues[5], result, 1e-3);
        Assert.assertEquals(6, generator.getNextIndex());

        for (int i = 6; i < referenceValues.length; i++) {
            result = generator.nextVector();
            Assert.assertArrayEquals(referenceValues[i], result, 1e-3);
            Assert.assertEquals(i + 1, generator.getNextIndex());
        }
    }

    @Test
    public void testSkip() {
        double[] result = generator.skipTo(5);
        Assert.assertArrayEquals(referenceValues[5], result, 1e-6);
        Assert.assertEquals(6, generator.getNextIndex());

        for (int i = 6; i < referenceValues.length; i++) {
            result = generator.nextVector();
            Assert.assertArrayEquals(referenceValues[i], result, 1e-6);
            Assert.assertEquals(i + 1, generator.getNextIndex());
        }
    }

[3318, 'src/test/java', 'org.apache.commons.math3.geometry.euclidean.threed', 'RotationTest', 'testComposeVectorOperator', 570, 586]

[3318, 'src/test/java', 'org.apache.commons.math3.geometry.euclidean.threed', 'RotationTest', 'testComposeInverseVectorOperator', 626, 642]

  @Test
  public void testComposeVectorOperator() throws MathIllegalArgumentException {

    Rotation r1 = new Rotation(new Vector3D(2, -3, 5), 1.7, RotationConvention.VECTOR_OPERATOR);
    Rotation r2 = new Rotation(new Vector3D(-1, 3, 2), 0.3, RotationConvention.VECTOR_OPERATOR);
    Rotation r3 = r2.compose(r1, RotationConvention.VECTOR_OPERATOR);

    for (double x = -0.9; x < 0.9; x += 0.2) {
      for (double y = -0.9; y < 0.9; y += 0.2) {
        for (double z = -0.9; z < 0.9; z += 0.2) {
          Vector3D u = new Vector3D(x, y, z);
          checkVector(r2.applyTo(r1.applyTo(u)), r3.applyTo(u));
        }
      }
    }

  }

  @Test
  public void testComposeInverseVectorOperator() throws MathIllegalArgumentException {

    Rotation r1 = new Rotation(new Vector3D(2, -3, 5), 1.7, RotationConvention.VECTOR_OPERATOR);
    Rotation r2 = new Rotation(new Vector3D(-1, 3, 2), 0.3, RotationConvention.VECTOR_OPERATOR);
    Rotation r3 = r2.composeInverse(r1, RotationConvention.VECTOR_OPERATOR);

    for (double x = -0.9; x < 0.9; x += 0.2) {
      for (double y = -0.9; y < 0.9; y += 0.2) {
        for (double z = -0.9; z < 0.9; z += 0.2) {
          Vector3D u = new Vector3D(x, y, z);
          checkVector(r2.applyInverseTo(r1.applyTo(u)), r3.applyTo(u));
        }
      }
    }

  }

[3344, 'src/test/java', 'org.apache.commons.math3.linear', 'RealVectorAbstractTest', 'testWalkInDefaultOrderPreservingVisitor2', 1416, 1464]

[3344, 'src/test/java', 'org.apache.commons.math3.linear', 'RealVectorAbstractTest', 'testWalkInOptimizedOrderPreservingVisitor2', 1539, 1587]

    @Test
    public void testWalkInDefaultOrderPreservingVisitor2() {
        final RealVector v = create(new double[5]);
        final RealVectorPreservingVisitor visitor;
        visitor = new RealVectorPreservingVisitor() {

            public void visit(int index, double value) {
                // Do nothing
            }

            public void start(int dimension, int start, int end) {
                // Do nothing
            }

            public double end() {
                return 0.0;
            }
        };
        try {
            v.walkInDefaultOrder(visitor, -1, 4);
            Assert.fail();
        } catch (OutOfRangeException e) {
            // Expected behavior
        }
        try {
            v.walkInDefaultOrder(visitor, 5, 4);
            Assert.fail();
        } catch (OutOfRangeException e) {
            // Expected behavior
        }
        try {
            v.walkInDefaultOrder(visitor, 0, -1);
            Assert.fail();
        } catch (OutOfRangeException e) {
            // Expected behavior
        }
        try {
            v.walkInDefaultOrder(visitor, 0, 5);
            Assert.fail();
        } catch (OutOfRangeException e) {
            // Expected behavior
        }
        try {
            v.walkInDefaultOrder(visitor, 4, 0);
            Assert.fail();
        } catch (NumberIsTooSmallException e) {
            // Expected behavior
        }
    }

    @Test
    public void testWalkInOptimizedOrderPreservingVisitor2() {
        final RealVector v = create(new double[5]);
        final RealVectorPreservingVisitor visitor;
        visitor = new RealVectorPreservingVisitor() {

            public void visit(int index, double value) {
                // Do nothing
            }

            public void start(int dimension, int start, int end) {
                // Do nothing
            }

            public double end() {
                return 0.0;
            }
        };
        try {
            v.walkInOptimizedOrder(visitor, -1, 4);
            Assert.fail();
        } catch (OutOfRangeException e) {
            // Expected behavior
        }
        try {
            v.walkInOptimizedOrder(visitor, 5, 4);
            Assert.fail();
        } catch (OutOfRangeException e) {
            // Expected behavior
        }
        try {
            v.walkInOptimizedOrder(visitor, 0, -1);
            Assert.fail();
        } catch (OutOfRangeException e) {
            // Expected behavior
        }
        try {
            v.walkInOptimizedOrder(visitor, 0, 5);
            Assert.fail();
        } catch (OutOfRangeException e) {
            // Expected behavior
        }
        try {
            v.walkInOptimizedOrder(visitor, 4, 0);
            Assert.fail();
        } catch (NumberIsTooSmallException e) {
            // Expected behavior
        }
    }

[3353, 'src/test/java', 'org.apache.commons.math3.analysis.solvers', 'MullerSolver2Test', 'testParameters', 125, 144]

[3353, 'src/test/java', 'org.apache.commons.math3.analysis.solvers', 'RiddersSolverTest', 'testParameters', 121, 140]

    @Test
    public void testParameters() {
        UnivariateFunction f = new Sin();
        UnivariateSolver solver = new MullerSolver2();

        try {
            // bad interval
            solver.solve(100, f, 1, -1);
            Assert.fail("Expecting NumberIsTooLargeException - bad interval");
        } catch (NumberIsTooLargeException ex) {
            // expected
        }
        try {
            // no bracketing
            solver.solve(100, f, 2, 3);
            Assert.fail("Expecting NoBracketingException - no bracketing");
        } catch (NoBracketingException ex) {
            // expected
        }
    }

    @Test
    public void testParameters() {
        UnivariateFunction f = new Sin();
        UnivariateSolver solver = new RiddersSolver();

        try {
            // bad interval
            solver.solve(100, f, 1, -1);
            Assert.fail("Expecting NumberIsTooLargeException - bad interval");
        } catch (NumberIsTooLargeException ex) {
            // expected
        }
        try {
            // no bracketing
            solver.solve(100, f, 2, 3);
            Assert.fail("Expecting NoBracketingException - no bracketing");
        } catch (NoBracketingException ex) {
            // expected
        }
    }

[3364, 'src/test/java', 'org.apache.commons.math3.analysis', 'FunctionUtilsTest', 'testSampleWrongBounds', 210, 213]

[3364, 'src/test/java', 'org.apache.commons.math3.analysis', 'FunctionUtilsTest', 'testSampleNullNumberOfPoints', 220, 223]

    @Test(expected = NumberIsTooLargeException.class)
    public void testSampleWrongBounds(){
        FunctionUtils.sample(new Sin(), FastMath.PI, 0.0, 10);
    }

    @Test(expected = NotStrictlyPositiveException.class)
    public void testSampleNullNumberOfPoints(){
        FunctionUtils.sample(new Sin(), 0.0, FastMath.PI, 0);
    }

[3397, 'src/test/java', 'org.apache.commons.math3.linear', 'SparseRealMatrixTest', 'testOperate', 252, 266]

[3397, 'src/test/java', 'org.apache.commons.math3.linear', 'SparseRealMatrixTest', 'testPremultiplyVector', 293, 307]

    @Test
    public void testOperate() {
        RealMatrix m = createSparseMatrix(id);
        assertClose("identity operate", testVector, m.operate(testVector),
                entryTolerance);
        assertClose("identity operate", testVector, m.operate(
                new ArrayRealVector(testVector)).toArray(), entryTolerance);
        m = createSparseMatrix(bigSingular);
        try {
            m.operate(testVector);
            Assert.fail("Expecting illegalArgumentException");
        } catch (MathIllegalArgumentException ex) {
            // ignored
        }
    }

    @Test
    public void testPremultiplyVector() {
        RealMatrix m = createSparseMatrix(testData);
        assertClose("premultiply", m.preMultiply(testVector), preMultTest,
            normTolerance);
        assertClose("premultiply", m.preMultiply(
            new ArrayRealVector(testVector).toArray()), preMultTest, normTolerance);
        m = createSparseMatrix(bigSingular);
        try {
            m.preMultiply(testVector);
            Assert.fail("expecting MathIllegalArgumentException");
        } catch (MathIllegalArgumentException ex) {
            // ignored
        }
    }

[3410, 'src/test/java', 'org.apache.commons.math3.linear', 'Array2DRowRealMatrixTest', 'testGetColumn', 833, 852]

[3410, 'src/test/java', 'org.apache.commons.math3.linear', 'BlockRealMatrixTest', 'testGetColumn', 972, 991]

    @Test
    public void testGetColumn() {
        RealMatrix m = new Array2DRowRealMatrix(subTestData);
        double[] mColumn1 = columnToArray(subColumn1);
        double[] mColumn3 = columnToArray(subColumn3);
        checkArrays(mColumn1, m.getColumn(1));
        checkArrays(mColumn3, m.getColumn(3));
        try {
            m.getColumn(-1);
            Assert.fail("Expecting OutOfRangeException");
        } catch (OutOfRangeException ex) {
            // expected
        }
        try {
            m.getColumn(4);
            Assert.fail("Expecting OutOfRangeException");
        } catch (OutOfRangeException ex) {
            // expected
        }
    }

    @Test
    public void testGetColumn() {
        RealMatrix m = new BlockRealMatrix(subTestData);
        double[] mColumn1 = columnToArray(subColumn1);
        double[] mColumn3 = columnToArray(subColumn3);
        checkArrays(mColumn1, m.getColumn(1));
        checkArrays(mColumn3, m.getColumn(3));
        try {
            m.getColumn(-1);
            Assert.fail("Expecting OutOfRangeException");
        } catch (OutOfRangeException ex) {
            // expected
        }
        try {
            m.getColumn(4);
            Assert.fail("Expecting OutOfRangeException");
        } catch (OutOfRangeException ex) {
            // expected
        }
    }

[3421, 'src/test/java', 'org.apache.commons.math3.analysis.interpolation', 'LoessInterpolatorTest', 'testNonStrictlyIncreasing1', 172, 175]

[3421, 'src/test/java', 'org.apache.commons.math3.analysis.interpolation', 'LoessInterpolatorTest', 'testNonStrictlyIncreasing2', 177, 180]

    @Test(expected=NonMonotonicSequenceException.class)
    public void testNonStrictlyIncreasing1() {
        new LoessInterpolator().smooth(new double[] {4,3,1,2}, new double[] {3,4,5,6});
    }

    @Test(expected=NonMonotonicSequenceException.class)
    public void testNonStrictlyIncreasing2() {
        new LoessInterpolator().smooth(new double[] {1,2,2,3}, new double[] {3,4,5,6});
    }

[3422, 'src/test/java', 'org.apache.commons.math3.analysis.function', 'LogisticTest', 'testParametricUsage5', 130, 134]

[3422, 'src/test/java', 'org.apache.commons.math3.analysis.function', 'LogisticTest', 'testParametricUsage6', 136, 140]

    @Test(expected=NotStrictlyPositiveException.class)
    public void testParametricUsage5() {
        final Logistic.Parametric g = new Logistic.Parametric();
        g.value(0, new double[] {1, 0, 1, 1, 0 ,0});
    }

    @Test(expected=NotStrictlyPositiveException.class)
    public void testParametricUsage6() {
        final Logistic.Parametric g = new Logistic.Parametric();
        g.gradient(0, new double[] {1, 0, 1, 1, 0 ,0});
    }

[3427, 'src/test/java', 'org.apache.commons.math3.fraction', 'BigFractionTest', 'testAbs', 291, 300]

[3427, 'src/test/java', 'org.apache.commons.math3.fraction', 'FractionTest', 'testAbs', 236, 245]

    @Test
    public void testAbs() {
        BigFraction a = new BigFraction(10, 21);
        BigFraction b = new BigFraction(-10, 21);
        BigFraction c = new BigFraction(10, -21);

        assertFraction(10, 21, a.abs());
        assertFraction(10, 21, b.abs());
        assertFraction(10, 21, c.abs());
    }

    @Test
    public void testAbs() {
        Fraction a = new Fraction(10, 21);
        Fraction b = new Fraction(-10, 21);
        Fraction c = new Fraction(10, -21);

        assertFraction(10, 21, a.abs());
        assertFraction(10, 21, b.abs());
        assertFraction(10, 21, c.abs());
    }

[3429, 'src/test/java', 'org.apache.commons.math3.complex', 'ComplexTest', 'testScalarAddInf', 147, 157]

[3429, 'src/test/java', 'org.apache.commons.math3.complex', 'ComplexTest', 'testScalarSubtractInf', 481, 490]

    @Test
    public void testScalarAddInf() {
        Complex x = new Complex(1, 1);
        double yDouble = Double.POSITIVE_INFINITY;

        Complex yComplex = new Complex(yDouble);
        Assert.assertEquals(x.add(yComplex), x.add(yDouble));

        x = new Complex(neginf, 0);
        Assert.assertEquals(x.add(yComplex), x.add(yDouble));
    }

    @Test
    public void testScalarSubtractInf() {
        Complex x = new Complex(1, 1);
        double yDouble = Double.POSITIVE_INFINITY;
        Complex yComplex = new Complex(yDouble);
        Assert.assertEquals(x.subtract(yComplex), x.subtract(yDouble));

        x = new Complex(neginf, 0);
        Assert.assertEquals(x.subtract(yComplex), x.subtract(yDouble));
    }

[3439, 'src/test/java', 'org.apache.commons.math3.analysis.interpolation', 'DividedDifferenceInterpolatorTest', 'testParameters', 114, 128]

[3439, 'src/test/java', 'org.apache.commons.math3.analysis.interpolation', 'NevilleInterpolatorTest', 'testParameters', 114, 128]

    @Test
    public void testParameters() {
        UnivariateInterpolator interpolator = new DividedDifferenceInterpolator();

        try {
            // bad abscissas array
            double x[] = { 1.0, 2.0, 2.0, 4.0 };
            double y[] = { 0.0, 4.0, 4.0, 2.5 };
            UnivariateFunction p = interpolator.interpolate(x, y);
            p.value(0.0);
            Assert.fail("Expecting NonMonotonicSequenceException - bad abscissas array");
        } catch (NonMonotonicSequenceException ex) {
            // expected
        }
    }

    @Test
    public void testParameters() {
        UnivariateInterpolator interpolator = new NevilleInterpolator();

        try {
            // bad abscissas array
            double x[] = { 1.0, 2.0, 2.0, 4.0 };
            double y[] = { 0.0, 4.0, 4.0, 2.5 };
            UnivariateFunction p = interpolator.interpolate(x, y);
            p.value(0.0);
            Assert.fail("Expecting NonMonotonicSequenceException - bad abscissas array");
        } catch (NonMonotonicSequenceException ex) {
            // expected
        }
    }

[3483, 'src/test/java', 'org.apache.commons.math3.stat.regression', 'MultipleLinearRegressionAbstractTest', 'canEstimateRegressionParameters', 42, 46]

[3483, 'src/test/java', 'org.apache.commons.math3.stat.regression', 'MultipleLinearRegressionAbstractTest', 'canEstimateResiduals', 48, 52]

    @Test
    public void canEstimateRegressionParameters(){
        double[] beta = regression.estimateRegressionParameters();
        Assert.assertEquals(getNumberOfRegressors(), beta.length);
    }

    @Test
    public void canEstimateResiduals(){
        double[] e = regression.estimateResiduals();
        Assert.assertEquals(getSampleSize(), e.length);
    }

[3495, 'src/test/java', 'org.apache.commons.math3.optimization.linear', 'SimplexSolverTest', 'testModelWithNoArtificialVars', 337, 350]

[3495, 'src/test/java', 'org.apache.commons.math3.optimization.linear', 'SimplexSolverTest', 'testSimplexSolver', 305, 319]

    @Test
    public void testModelWithNoArtificialVars() {
        LinearObjectiveFunction f = new LinearObjectiveFunction(new double[] { 15, 10 }, 0);
        Collection<LinearConstraint> constraints = new ArrayList<LinearConstraint>();
        constraints.add(new LinearConstraint(new double[] { 1, 0 }, Relationship.LEQ, 2));
        constraints.add(new LinearConstraint(new double[] { 0, 1 }, Relationship.LEQ, 3));
        constraints.add(new LinearConstraint(new double[] { 1, 1 }, Relationship.LEQ, 4));

        SimplexSolver solver = new SimplexSolver();
        PointValuePair solution = solver.optimize(f, constraints, GoalType.MAXIMIZE, false);
        Assert.assertEquals(2.0, solution.getPoint()[0], 0.0);
        Assert.assertEquals(2.0, solution.getPoint()[1], 0.0);
        Assert.assertEquals(50.0, solution.getValue(), 0.0);
    }

    @Test
    public void testSimplexSolver() {
        LinearObjectiveFunction f =
            new LinearObjectiveFunction(new double[] { 15, 10 }, 7);
        Collection<LinearConstraint> constraints = new ArrayList<LinearConstraint>();
        constraints.add(new LinearConstraint(new double[] { 1, 0 }, Relationship.LEQ, 2));
        constraints.add(new LinearConstraint(new double[] { 0, 1 }, Relationship.LEQ, 3));
        constraints.add(new LinearConstraint(new double[] { 1, 1 }, Relationship.EQ, 4));

        SimplexSolver solver = new SimplexSolver();
        PointValuePair solution = solver.optimize(f, constraints, GoalType.MAXIMIZE, false);
        Assert.assertEquals(2.0, solution.getPoint()[0], 0.0);
        Assert.assertEquals(2.0, solution.getPoint()[1], 0.0);
        Assert.assertEquals(57.0, solution.getValue(), 0.0);
    }

[3505, 'src/test/java', 'org.apache.commons.math3.linear', 'FieldLUSolverTest', 'testSolveDimensionErrors', 80, 98]

[3505, 'src/test/java', 'org.apache.commons.math3.linear', 'FieldLUSolverTest', 'testSolveSingularityErrors', 101, 119]

    @Test
    public void testSolveDimensionErrors() {
        FieldDecompositionSolver<Fraction> solver;
        solver = new FieldLUDecomposition<Fraction>(createFractionMatrix(testData))
            .getSolver();
        FieldMatrix<Fraction> b = createFractionMatrix(new int[2][2]);
        try {
            solver.solve(b);
            Assert.fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException iae) {
            // expected behavior
        }
        try {
            solver.solve(b.getColumnVector(0));
            Assert.fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException iae) {
            // expected behavior
        }
    }

    @Test
    public void testSolveSingularityErrors() {
        FieldDecompositionSolver<Fraction> solver;
        solver = new FieldLUDecomposition<Fraction>(createFractionMatrix(singular))
            .getSolver();
        FieldMatrix<Fraction> b = createFractionMatrix(new int[2][2]);
        try {
            solver.solve(b);
            Assert.fail("an exception should have been thrown");
        } catch (SingularMatrixException ime) {
            // expected behavior
        }
        try {
            solver.solve(b.getColumnVector(0));
            Assert.fail("an exception should have been thrown");
        } catch (SingularMatrixException ime) {
            // expected behavior
        }
    }

[3507, 'src/test/java', 'org.apache.commons.math3.fraction', 'BigFractionTest', 'testDoubleValueForLargeNumeratorAndDenominator', 218, 227]

[3507, 'src/test/java', 'org.apache.commons.math3.fraction', 'BigFractionTest', 'testFloatValueForLargeNumeratorAndDenominator', 230, 239]

    @Test
    public void testDoubleValueForLargeNumeratorAndDenominator() {
        final BigInteger pow400 = BigInteger.TEN.pow(400);
        final BigInteger pow401 = BigInteger.TEN.pow(401);
        final BigInteger two = new BigInteger("2");
        final BigFraction large = new BigFraction(pow401.add(BigInteger.ONE),
                                                  pow400.multiply(two));

        Assert.assertEquals(5, large.doubleValue(), 1e-15);
    }

    @Test
    public void testFloatValueForLargeNumeratorAndDenominator() {
        final BigInteger pow400 = BigInteger.TEN.pow(400);
        final BigInteger pow401 = BigInteger.TEN.pow(401);
        final BigInteger two = new BigInteger("2");
        final BigFraction large = new BigFraction(pow401.add(BigInteger.ONE),
                                                  pow400.multiply(two));

        Assert.assertEquals(5, large.floatValue(), 1e-15);
    }

[3508, 'src/test/java', 'org.apache.commons.math3.stat.descriptive.rank', 'PSquarePercentileTest', 'testMarkerHeightWithLowerIndex', 298, 306]

[3508, 'src/test/java', 'org.apache.commons.math3.stat.descriptive.rank', 'PSquarePercentileTest', 'testMarkerHeightWithHigherIndex', 308, 316]

    @Test(expected = OutOfRangeException.class)
    public void testMarkerHeightWithLowerIndex() {
        PSquareMarkers mThat =
                PSquarePercentile.newMarkers(
                        Arrays.asList(new Double[] { 95.1772, 95.1567, 95.1937,
                                95.1959, 95.1442, 95.0610, 95.1591, 95.1195,
                                95.1772, 95.0925, 95.1990, 95.1682 }), 0.50);
        mThat.height(0);
    }

    @Test(expected = OutOfRangeException.class)
    public void testMarkerHeightWithHigherIndex() {
        PSquareMarkers mThat =
                PSquarePercentile.newMarkers(
                        Arrays.asList(new Double[] { 95.1772, 95.1567, 95.1937,
                                95.1959, 95.1442, 95.0610, 95.1591, 95.1195,
                                95.1772, 95.0925, 95.1990, 95.1682 }), 0.50);
        mThat.height(6);
    }

[3509, 'src/test/java', 'org.apache.commons.math3.geometry.spherical.oned', 'ArcsSetTest', 'testFullEqualEndPoints', 119, 131]

[3509, 'src/test/java', 'org.apache.commons.math3.geometry.spherical.oned', 'ArcsSetTest', 'testFullCircle', 133, 145]

    @Test
    public void testFullEqualEndPoints() {
        ArcsSet set = new ArcsSet(1.0, 1.0, 1.0e-10);
        Assert.assertEquals(1.0e-10, set.getTolerance(), 1.0e-20);
        Assert.assertEquals(Region.Location.INSIDE, set.checkPoint(new S1Point(9.0)));
        for (double alpha = -20.0; alpha <= 20.0; alpha += 0.1) {
            Assert.assertEquals(Region.Location.INSIDE, set.checkPoint(new S1Point(alpha)));
        }
        Assert.assertEquals(1, set.asList().size());
        Assert.assertEquals(0.0, set.asList().get(0).getInf(), 1.0e-10);
        Assert.assertEquals(2 * FastMath.PI, set.asList().get(0).getSup(), 1.0e-10);
        Assert.assertEquals(2 * FastMath.PI, set.getSize(), 1.0e-10);
    }

    @Test
    public void testFullCircle() {
        ArcsSet set = new ArcsSet(1.0e-10);
        Assert.assertEquals(1.0e-10, set.getTolerance(), 1.0e-20);
        Assert.assertEquals(Region.Location.INSIDE, set.checkPoint(new S1Point(9.0)));
        for (double alpha = -20.0; alpha <= 20.0; alpha += 0.1) {
            Assert.assertEquals(Region.Location.INSIDE, set.checkPoint(new S1Point(alpha)));
        }
        Assert.assertEquals(1, set.asList().size());
        Assert.assertEquals(0.0, set.asList().get(0).getInf(), 1.0e-10);
        Assert.assertEquals(2 * FastMath.PI, set.asList().get(0).getSup(), 1.0e-10);
        Assert.assertEquals(2 * FastMath.PI, set.getSize(), 1.0e-10);
    }

[3551, 'src/test/java', 'org.apache.commons.math3.random', 'GaussianRandomGeneratorTest', 'testMeanAndStandardDeviation', 27, 38]

[3551, 'src/test/java', 'org.apache.commons.math3.random', 'UniformRandomGeneratorTest', 'testMeanAndStandardDeviation', 27, 38]

    @Test
    public void testMeanAndStandardDeviation() {
        RandomGenerator rg = new JDKRandomGenerator();
        rg.setSeed(17399225432l);
        GaussianRandomGenerator generator = new GaussianRandomGenerator(rg);
        double[] sample = new double[10000];
        for (int i = 0; i < sample.length; ++i) {
            sample[i] = generator.nextNormalizedDouble();
        }
        Assert.assertEquals(0.0, StatUtils.mean(sample), 0.012);
        Assert.assertEquals(1.0, StatUtils.variance(sample), 0.01);
    }

    @Test
    public void testMeanAndStandardDeviation() {
        RandomGenerator rg = new JDKRandomGenerator();
        rg.setSeed(17399225432l);
        UniformRandomGenerator generator = new UniformRandomGenerator(rg);
        double[] sample = new double[10000];
        for (int i = 0; i < sample.length; ++i) {
            sample[i] = generator.nextNormalizedDouble();
        }
        Assert.assertEquals(0.0, StatUtils.mean(sample), 0.07);
        Assert.assertEquals(1.0, StatUtils.variance(sample), 0.02);
    }

[3553, 'src/test/java', 'org.apache.commons.math3.linear', 'BlockFieldMatrixTest', 'testGetRow', 1001, 1018]

[3553, 'src/test/java', 'org.apache.commons.math3.linear', 'FieldMatrixImplTest', 'testGetRow', 761, 778]

    @Test
    public void testGetRow() {
        FieldMatrix<Fraction> m = new BlockFieldMatrix<Fraction>(subTestData);
        checkArrays(subRow0[0], m.getRow(0));
        checkArrays(subRow3[0], m.getRow(3));
        try {
            m.getRow(-1);
            Assert.fail("Expecting OutOfRangeException");
        } catch (OutOfRangeException ex) {
            // expected
        }
        try {
            m.getRow(4);
            Assert.fail("Expecting OutOfRangeException");
        } catch (OutOfRangeException ex) {
            // expected
        }
    }

    @Test
    public void testGetRow() {
        FieldMatrix<Fraction> m = new Array2DRowFieldMatrix<Fraction>(subTestData);
        checkArrays(subRow0[0], m.getRow(0));
        checkArrays(subRow3[0], m.getRow(3));
        try {
            m.getRow(-1);
            Assert.fail("Expecting OutOfRangeException");
        } catch (OutOfRangeException ex) {
            // expected
        }
        try {
            m.getRow(4);
            Assert.fail("Expecting OutOfRangeException");
        } catch (OutOfRangeException ex) {
            // expected
        }
    }

[3559, 'src/test/java', 'org.apache.commons.math3.analysis.integration.gauss', 'LegendreHighPrecisionTest', 'testCos', 33, 41]

[3559, 'src/test/java', 'org.apache.commons.math3.analysis.integration.gauss', 'LegendreTest', 'testCos', 33, 41]

    @Test
    public void testCos() {
        final UnivariateFunction cos = new Cos();

        final GaussIntegrator integrator = factory.legendreHighPrecision(7, 0, Math.PI / 2);
        final double s = integrator.integrate(cos);
        // System.out.println("s=" + s + " e=" + 1);
        Assert.assertEquals(1, s, Math.ulp(1d));
    }

    @Test
    public void testCos() {
        final UnivariateFunction cos = new Cos();

        final GaussIntegrator integrator = factory.legendre(7, 0, Math.PI / 2);
        final double s = integrator.integrate(cos);
        // System.out.println("s=" + s + " e=" + 1);
        Assert.assertEquals(1, s, Math.ulp(1d));
    }

[3573, 'src/test/java', 'org.apache.commons.math3.linear', 'BlockFieldMatrixTest', 'testSetRowVector', 890, 909]

[3573, 'src/test/java', 'org.apache.commons.math3.linear', 'FieldMatrixImplTest', 'testSetRowVector', 690, 709]

    @Test
    public void testSetRowVector() {
        FieldMatrix<Fraction> m = new BlockFieldMatrix<Fraction>(subTestData);
        FieldVector<Fraction> mRow3 = new ArrayFieldVector<Fraction>(subRow3[0]);
        Assert.assertNotSame(mRow3, m.getRowMatrix(0));
        m.setRowVector(0, mRow3);
        Assert.assertEquals(mRow3, m.getRowVector(0));
        try {
            m.setRowVector(-1, mRow3);
            Assert.fail("Expecting OutOfRangeException");
        } catch (OutOfRangeException ex) {
            // expected
        }
        try {
            m.setRowVector(0, new ArrayFieldVector<Fraction>(FractionField.getInstance(), 5));
            Assert.fail("Expecting MatrixDimensionMismatchException");
        } catch (MatrixDimensionMismatchException ex) {
            // expected
        }
    }

    @Test
    public void testSetRowVector() {
        FieldMatrix<Fraction> m = new Array2DRowFieldMatrix<Fraction>(subTestData);
        FieldVector<Fraction> mRow3 = new ArrayFieldVector<Fraction>(subRow3[0]);
        Assert.assertNotSame(mRow3, m.getRowMatrix(0));
        m.setRowVector(0, mRow3);
        Assert.assertEquals(mRow3, m.getRowVector(0));
        try {
            m.setRowVector(-1, mRow3);
            Assert.fail("Expecting OutOfRangeException");
        } catch (OutOfRangeException ex) {
            // expected
        }
        try {
            m.setRowVector(0, new ArrayFieldVector<Fraction>(FractionField.getInstance(), 5));
            Assert.fail("Expecting MatrixDimensionMismatchException");
        } catch (MatrixDimensionMismatchException ex) {
            // expected
        }
    }

[3582, 'src/test/java', 'org.apache.commons.math3.linear', 'RealVectorAbstractTest', 'testGetSubVectorInvalidIndex1', 371, 375]

[3582, 'src/test/java', 'org.apache.commons.math3.linear', 'RealVectorAbstractTest', 'testGetSubVectorInvalidIndex4', 389, 393]

    @Test(expected = OutOfRangeException.class)
    public void testGetSubVectorInvalidIndex1() {
        final int n = 10;
        create(new double[n]).getSubVector(-1, 2);
    }

    @Test(expected = NotPositiveException.class)
    public void testGetSubVectorInvalidIndex4() {
        final int n = 10;
        create(new double[n]).getSubVector(3, -2);
    }

[3583, 'src/test/java', 'org.apache.commons.math3.linear', 'HessenbergTransformerTest', 'testAEqualPHPt', 61, 66]

[3583, 'src/test/java', 'org.apache.commons.math3.linear', 'SchurTransformerTest', 'testAEqualPTPt', 62, 67]

    @Test
    public void testAEqualPHPt() {
        checkAEqualPHPt(MatrixUtils.createRealMatrix(testSquare5));
        checkAEqualPHPt(MatrixUtils.createRealMatrix(testSquare3));
        checkAEqualPHPt(MatrixUtils.createRealMatrix(testRandom));
   }

    @Test
    public void testAEqualPTPt() {
        checkAEqualPTPt(MatrixUtils.createRealMatrix(testSquare5));
        checkAEqualPTPt(MatrixUtils.createRealMatrix(testSquare3));
        checkAEqualPTPt(MatrixUtils.createRealMatrix(testRandom));
   }

[3599, 'src/test/java', 'org.apache.commons.math3.geometry.euclidean.oned', 'Vector1DFormatAbstractTest', 'testSimpleWithDecimals', 52, 60]

[3599, 'src/test/java', 'org.apache.commons.math3.geometry.euclidean.oned', 'Vector1DFormatAbstractTest', 'testSimpleWithDecimalsTrunc', 62, 70]

    @Test
    public void testSimpleWithDecimals() {
        Vector1D c = new Vector1D(1.23);
        String expected =
            "{1"    + getDecimalCharacter() +
            "23}";
        String actual = vector1DFormat.format(c);
        Assert.assertEquals(expected, actual);
    }

    @Test
    public void testSimpleWithDecimalsTrunc() {
        Vector1D c = new Vector1D(1.232323232323);
        String expected =
            "{1"    + getDecimalCharacter() +
            "2323232323}";
        String actual = vector1DFormat.format(c);
        Assert.assertEquals(expected, actual);
    }

[3616, 'src/test/java', 'org.apache.commons.math3.random', 'RandomDataGeneratorTest', 'testNextChiSquare', 1036, 1046]

[3616, 'src/test/java', 'org.apache.commons.math3.random', 'RandomDataGeneratorTest', 'testNextT', 1086, 1096]

    @Test
    public void testNextChiSquare() {
        double[] quartiles = TestUtils.getDistributionQuartiles(new ChiSquaredDistribution(12));
        long[] counts = new long[4];
        randomData.reSeed(1000);
        for (int i = 0; i < 1000; i++) {
            double value = randomData.nextChiSquare(12);
            TestUtils.updateCounts(value, counts, quartiles);
        }
        TestUtils.assertChiSquareAccept(expected, counts, 0.001);
    }

    @Test
    public void testNextT() {
        double[] quartiles = TestUtils.getDistributionQuartiles(new TDistribution(10));
        long[] counts = new long[4];
        randomData.reSeed(1000);
        for (int i = 0; i < 1000; i++) {
            double value = randomData.nextT(10);
            TestUtils.updateCounts(value, counts, quartiles);
        }
        TestUtils.assertChiSquareAccept(expected, counts, 0.001);
    }

[3618, 'src/test/java', 'org.apache.commons.math3.linear', 'BlockFieldMatrixTest', 'testPremultiply', 460, 480]

[3618, 'src/test/java', 'org.apache.commons.math3.linear', 'FieldMatrixImplTest', 'testPremultiply', 323, 343]

    @Test
    public void testPremultiply() {
        FieldMatrix<Fraction> m3 = new BlockFieldMatrix<Fraction>(d3);
        FieldMatrix<Fraction> m4 = new BlockFieldMatrix<Fraction>(d4);
        FieldMatrix<Fraction> m5 = new BlockFieldMatrix<Fraction>(d5);
        TestUtils.assertEquals(m4.preMultiply(m3), m5);

        BlockFieldMatrix<Fraction> m = new BlockFieldMatrix<Fraction>(testData);
        BlockFieldMatrix<Fraction> mInv = new BlockFieldMatrix<Fraction>(testDataInv);
        BlockFieldMatrix<Fraction> identity = new BlockFieldMatrix<Fraction>(id);
        TestUtils.assertEquals(m.preMultiply(mInv), identity);
        TestUtils.assertEquals(mInv.preMultiply(m), identity);
        TestUtils.assertEquals(m.preMultiply(identity), m);
        TestUtils.assertEquals(identity.preMultiply(mInv), mInv);
        try {
            m.preMultiply(new BlockFieldMatrix<Fraction>(bigSingular));
            Assert.fail("Expecting illegalArgumentException");
        } catch (MathIllegalArgumentException ex) {
            // ignored
        }
    }

    @Test
    public void testPremultiply() {
        FieldMatrix<Fraction> m3 = new Array2DRowFieldMatrix<Fraction>(d3);
        FieldMatrix<Fraction> m4 = new Array2DRowFieldMatrix<Fraction>(d4);
        FieldMatrix<Fraction> m5 = new Array2DRowFieldMatrix<Fraction>(d5);
        TestUtils.assertEquals(m4.preMultiply(m3), m5);

        Array2DRowFieldMatrix<Fraction> m = new Array2DRowFieldMatrix<Fraction>(testData);
        Array2DRowFieldMatrix<Fraction> mInv = new Array2DRowFieldMatrix<Fraction>(testDataInv);
        Array2DRowFieldMatrix<Fraction> identity = new Array2DRowFieldMatrix<Fraction>(id);
        TestUtils.assertEquals(m.preMultiply(mInv), identity);
        TestUtils.assertEquals(mInv.preMultiply(m), identity);
        TestUtils.assertEquals(m.preMultiply(identity), m);
        TestUtils.assertEquals(identity.preMultiply(mInv), mInv);
        try {
            m.preMultiply(new Array2DRowFieldMatrix<Fraction>(bigSingular));
            Assert.fail("Expecting illegalArgumentException");
        } catch (MathIllegalArgumentException ex) {
            // ignored
        }
    }

[3663, 'src/test/java', 'org.apache.commons.math3.linear', 'Array2DRowRealMatrixTest', 'testGetRow', 794, 811]

[3663, 'src/test/java', 'org.apache.commons.math3.linear', 'BlockRealMatrixTest', 'testGetRow', 913, 930]

    @Test
    public void testGetRow() {
        RealMatrix m = new Array2DRowRealMatrix(subTestData);
        checkArrays(subRow0[0], m.getRow(0));
        checkArrays(subRow3[0], m.getRow(3));
        try {
            m.getRow(-1);
            Assert.fail("Expecting OutOfRangeException");
        } catch (OutOfRangeException ex) {
            // expected
        }
        try {
            m.getRow(4);
            Assert.fail("Expecting OutOfRangeException");
        } catch (OutOfRangeException ex) {
            // expected
        }
    }

    @Test
    public void testGetRow() {
        RealMatrix m = new BlockRealMatrix(subTestData);
        checkArrays(subRow0[0], m.getRow(0));
        checkArrays(subRow3[0], m.getRow(3));
        try {
            m.getRow(-1);
            Assert.fail("Expecting OutOfRangeException");
        } catch (OutOfRangeException ex) {
            // expected
        }
        try {
            m.getRow(4);
            Assert.fail("Expecting OutOfRangeException");
        } catch (OutOfRangeException ex) {
            // expected
        }
    }

[3665, 'src/test/java', 'org.apache.commons.math3.fraction', 'BigFractionTest', 'testSerial', 628, 638]

[3665, 'src/test/java', 'org.apache.commons.math3.fraction', 'FractionTest', 'testSerial', 617, 627]

    @Test
    public void testSerial() throws FractionConversionException {
        BigFraction[] fractions = {
            new BigFraction(3, 4), BigFraction.ONE, BigFraction.ZERO,
            new BigFraction(17), new BigFraction(FastMath.PI, 1000),
            new BigFraction(-5, 2)
        };
        for (BigFraction fraction : fractions) {
            Assert.assertEquals(fraction, TestUtils.serializeAndRecover(fraction));
        }
    }

    @Test
    public void testSerial() throws FractionConversionException {
        Fraction[] fractions = {
            new Fraction(3, 4), Fraction.ONE, Fraction.ZERO,
            new Fraction(17), new Fraction(FastMath.PI, 1000),
            new Fraction(-5, 2)
        };
        for (Fraction fraction : fractions) {
            Assert.assertEquals(fraction, TestUtils.serializeAndRecover(fraction));
        }
    }

