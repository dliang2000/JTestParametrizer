[21, 'src/test/java', 'org.apache.commons.math3.optimization.direct', 'SimplexOptimizerMultiDirectionalTest', 'testMinimize2', 45, 58]

[21, 'src/test/java', 'org.apache.commons.math3.optimization.direct', 'SimplexOptimizerNelderMeadTest', 'testMinimize2', 50, 63]

[21, 'src/test/java', 'org.apache.commons.math3.optimization.direct', 'SimplexOptimizerNelderMeadTest', 'testMaximize2', 80, 93]

    @Test
    public void testMinimize2() {
        SimplexOptimizer optimizer = new SimplexOptimizer(1e-11, 1e-30);
        optimizer.setSimplex(new MultiDirectionalSimplex(new double[] { 0.2, 0.2 }));
        final FourExtrema fourExtrema = new FourExtrema();

        final PointValuePair optimum
            =  optimizer.optimize(200, fourExtrema, GoalType.MINIMIZE, new double[] { 1, 0 });
        Assert.assertEquals(fourExtrema.xP, optimum.getPoint()[0], 2e-8);
        Assert.assertEquals(fourExtrema.yM, optimum.getPoint()[1], 3e-6);
        Assert.assertEquals(fourExtrema.valueXpYm, optimum.getValue(), 2e-12);
        Assert.assertTrue(optimizer.getEvaluations() > 120);
        Assert.assertTrue(optimizer.getEvaluations() < 150);
    }

    @Test
    public void testMinimize2() {
        SimplexOptimizer optimizer = new SimplexOptimizer(1e-10, 1e-30);
        optimizer.setSimplex(new NelderMeadSimplex(new double[] { 0.2, 0.2 }));
        final FourExtrema fourExtrema = new FourExtrema();

        final PointValuePair optimum
            = optimizer.optimize(100, fourExtrema, GoalType.MINIMIZE, new double[] { 1, 0 });
        Assert.assertEquals(fourExtrema.xP, optimum.getPoint()[0], 5e-6);
        Assert.assertEquals(fourExtrema.yM, optimum.getPoint()[1], 6e-6);
        Assert.assertEquals(fourExtrema.valueXpYm, optimum.getValue(), 1e-11);
        Assert.assertTrue(optimizer.getEvaluations() > 60);
        Assert.assertTrue(optimizer.getEvaluations() < 90);
    }

    @Test
    public void testMaximize2() {
        SimplexOptimizer optimizer = new SimplexOptimizer(1e-10, 1e-30);
        optimizer.setSimplex(new NelderMeadSimplex(new double[] { 0.2, 0.2 }));
        final FourExtrema fourExtrema = new FourExtrema();

        final PointValuePair optimum
            = optimizer.optimize(100, fourExtrema, GoalType.MAXIMIZE, new double[] { 1, 0 });
        Assert.assertEquals(fourExtrema.xP, optimum.getPoint()[0], 4e-6);
        Assert.assertEquals(fourExtrema.yP, optimum.getPoint()[1], 5e-6);
        Assert.assertEquals(fourExtrema.valueXpYp, optimum.getValue(), 7e-12);
        Assert.assertTrue(optimizer.getEvaluations() > 60);
        Assert.assertTrue(optimizer.getEvaluations() < 90);
    }

[49, 'src/test/java', 'org.apache.commons.math3.ode.nonstiff', 'ClassicalRungeKuttaIntegratorTest', 'testBackward', 224, 242]

[49, 'src/test/java', 'org.apache.commons.math3.ode.nonstiff', 'GillIntegratorTest', 'testBackward', 144, 162]

[49, 'src/test/java', 'org.apache.commons.math3.ode.nonstiff', 'LutherIntegratorTest', 'testBackward', 224, 242]

[49, 'src/test/java', 'org.apache.commons.math3.ode.nonstiff', 'MidpointIntegratorTest', 'testBackward', 146, 164]

[49, 'src/test/java', 'org.apache.commons.math3.ode.nonstiff', 'ThreeEighthesIntegratorTest', 'testBackward', 144, 162]

  @Test
  public void testBackward()
      throws DimensionMismatchException, NumberIsTooSmallException,
             MaxCountExceededException, NoBracketingException {

    TestProblem5 pb = new TestProblem5();
    double step = FastMath.abs(pb.getFinalTime() - pb.getInitialTime()) * 0.001;

    FirstOrderIntegrator integ = new ClassicalRungeKuttaIntegrator(step);
    TestProblemHandler handler = new TestProblemHandler(pb, integ);
    integ.addStepHandler(handler);
    integ.integrate(pb, pb.getInitialTime(), pb.getInitialState(),
                    pb.getFinalTime(), new double[pb.getDimension()]);

    Assert.assertTrue(handler.getLastError() < 5.0e-10);
    Assert.assertTrue(handler.getMaximalValueError() < 7.0e-10);
    Assert.assertEquals(0, handler.getMaximalTimeError(), 1.0e-12);
    Assert.assertEquals("classical Runge-Kutta", integ.getName());
  }

  @Test
  public void testBackward()
      throws DimensionMismatchException, NumberIsTooSmallException,
             MaxCountExceededException, NoBracketingException {

      TestProblem5 pb = new TestProblem5();
      double step = FastMath.abs(pb.getFinalTime() - pb.getInitialTime()) * 0.001;

      FirstOrderIntegrator integ = new GillIntegrator(step);
      TestProblemHandler handler = new TestProblemHandler(pb, integ);
      integ.addStepHandler(handler);
      integ.integrate(pb, pb.getInitialTime(), pb.getInitialState(),
                      pb.getFinalTime(), new double[pb.getDimension()]);

      Assert.assertTrue(handler.getLastError() < 5.0e-10);
      Assert.assertTrue(handler.getMaximalValueError() < 7.0e-10);
      Assert.assertEquals(0, handler.getMaximalTimeError(), 1.0e-12);
      Assert.assertEquals("Gill", integ.getName());
  }

    @Test
    public void testBackward()
            throws DimensionMismatchException, NumberIsTooSmallException,
            MaxCountExceededException, NoBracketingException {

        TestProblem5 pb = new TestProblem5();
        double step = FastMath.abs(pb.getFinalTime() - pb.getInitialTime()) * 0.001;

        FirstOrderIntegrator integ = new LutherIntegrator(step);
        TestProblemHandler handler = new TestProblemHandler(pb, integ);
        integ.addStepHandler(handler);
        integ.integrate(pb, pb.getInitialTime(), pb.getInitialState(),
                        pb.getFinalTime(), new double[pb.getDimension()]);

        Assert.assertTrue(handler.getLastError() < 3.0e-13);
        Assert.assertTrue(handler.getMaximalValueError() < 5.0e-13);
        Assert.assertEquals(0, handler.getMaximalTimeError(), 1.0e-12);
        Assert.assertEquals("Luther", integ.getName());
    }

  @Test
  public void testBackward()
      throws DimensionMismatchException, NumberIsTooSmallException,
             MaxCountExceededException, NoBracketingException {

      TestProblem5 pb = new TestProblem5();
      double step = FastMath.abs(pb.getFinalTime() - pb.getInitialTime()) * 0.001;

      FirstOrderIntegrator integ = new MidpointIntegrator(step);
      TestProblemHandler handler = new TestProblemHandler(pb, integ);
      integ.addStepHandler(handler);
      integ.integrate(pb, pb.getInitialTime(), pb.getInitialState(),
                      pb.getFinalTime(), new double[pb.getDimension()]);

      Assert.assertTrue(handler.getLastError() < 6.0e-4);
      Assert.assertTrue(handler.getMaximalValueError() < 6.0e-4);
      Assert.assertEquals(0, handler.getMaximalTimeError(), 1.0e-12);
      Assert.assertEquals("midpoint", integ.getName());
  }

  @Test
  public void testBackward()
      throws DimensionMismatchException, NumberIsTooSmallException,
             MaxCountExceededException, NoBracketingException {

      TestProblem5 pb = new TestProblem5();
      double step = FastMath.abs(pb.getFinalTime() - pb.getInitialTime()) * 0.001;

      FirstOrderIntegrator integ = new ThreeEighthesIntegrator(step);
      TestProblemHandler handler = new TestProblemHandler(pb, integ);
      integ.addStepHandler(handler);
      integ.integrate(pb, pb.getInitialTime(), pb.getInitialState(),
                      pb.getFinalTime(), new double[pb.getDimension()]);

      Assert.assertTrue(handler.getLastError() < 5.0e-10);
      Assert.assertTrue(handler.getMaximalValueError() < 7.0e-10);
      Assert.assertEquals(0, handler.getMaximalTimeError(), 1.0e-12);
      Assert.assertEquals("3/8", integ.getName());
  }

[85, 'src/test/java', 'org.apache.commons.math3.special', 'BetaTest', 'testLogGammaSum', 317, 331]

[85, 'src/test/java', 'org.apache.commons.math3.special', 'BetaTest', 'testLogGammaMinusLogGammaSum', 503, 517]

[85, 'src/test/java', 'org.apache.commons.math3.special', 'BetaTest', 'testSumDeltaMinusDeltaSum', 676, 691]

    @Test
    public void testLogGammaSum() {
        final int ulps = 2;
        for (int i = 0; i < LOG_GAMMA_SUM_REF.length; i++) {
            final double[] ref = LOG_GAMMA_SUM_REF[i];
            final double a = ref[0];
            final double b = ref[1];
            final double expected = ref[2];
            final double actual = logGammaSum(a, b);
            final double tol = ulps * FastMath.ulp(expected);
            final StringBuilder builder = new StringBuilder();
            builder.append(a).append(", ").append(b);
            Assert.assertEquals(builder.toString(), expected, actual, tol);
        }
    }

    @Test
    public void testLogGammaMinusLogGammaSum() {
        final int ulps = 4;
        for (int i = 0; i < LOG_GAMMA_MINUS_LOG_GAMMA_SUM_REF.length; i++) {
            final double[] ref = LOG_GAMMA_MINUS_LOG_GAMMA_SUM_REF[i];
            final double a = ref[0];
            final double b = ref[1];
            final double expected = ref[2];
            final double actual = logGammaMinusLogGammaSum(a, b);
            final double tol = ulps * FastMath.ulp(expected);
            final StringBuilder builder = new StringBuilder();
            builder.append(a).append(", ").append(b);
            Assert.assertEquals(builder.toString(), expected, actual, tol);
        }
    }

    @Test
    public void testSumDeltaMinusDeltaSum() {

        final int ulps = 3;
        for (int i = 0; i < SUM_DELTA_MINUS_DELTA_SUM_REF.length; i++) {
            final double[] ref = SUM_DELTA_MINUS_DELTA_SUM_REF[i];
            final double a = ref[0];
            final double b = ref[1];
            final double expected = ref[2];
            final double actual = sumDeltaMinusDeltaSum(a, b);
            final double tol = ulps * FastMath.ulp(expected);
            final StringBuilder builder = new StringBuilder();
            builder.append(a).append(", ").append(b);
            Assert.assertEquals(builder.toString(), expected, actual, tol);
        }
    }

[160, 'src/test/java', 'org.apache.commons.math3.analysis.solvers', 'MullerSolver2Test', 'testSinFunction', 45, 62]

[160, 'src/test/java', 'org.apache.commons.math3.analysis.solvers', 'MullerSolverTest', 'testSinFunction', 45, 62]

[160, 'src/test/java', 'org.apache.commons.math3.analysis.solvers', 'RiddersSolverTest', 'testSinFunction', 43, 60]

    @Test
    public void testSinFunction() {
        UnivariateFunction f = new Sin();
        UnivariateSolver solver = new MullerSolver2();
        double min, max, expected, result, tolerance;

        min = 3.0; max = 4.0; expected = FastMath.PI;
        tolerance = FastMath.max(solver.getAbsoluteAccuracy(),
                    FastMath.abs(expected * solver.getRelativeAccuracy()));
        result = solver.solve(100, f, min, max);
        Assert.assertEquals(expected, result, tolerance);

        min = -1.0; max = 1.5; expected = 0.0;
        tolerance = FastMath.max(solver.getAbsoluteAccuracy(),
                    FastMath.abs(expected * solver.getRelativeAccuracy()));
        result = solver.solve(100, f, min, max);
        Assert.assertEquals(expected, result, tolerance);
    }

    @Test
    public void testSinFunction() {
        UnivariateFunction f = new Sin();
        UnivariateSolver solver = new MullerSolver();
        double min, max, expected, result, tolerance;

        min = 3.0; max = 4.0; expected = FastMath.PI;
        tolerance = FastMath.max(solver.getAbsoluteAccuracy(),
                    FastMath.abs(expected * solver.getRelativeAccuracy()));
        result = solver.solve(100, f, min, max);
        Assert.assertEquals(expected, result, tolerance);

        min = -1.0; max = 1.5; expected = 0.0;
        tolerance = FastMath.max(solver.getAbsoluteAccuracy(),
                    FastMath.abs(expected * solver.getRelativeAccuracy()));
        result = solver.solve(100, f, min, max);
        Assert.assertEquals(expected, result, tolerance);
    }

    @Test
    public void testSinFunction() {
        UnivariateFunction f = new Sin();
        UnivariateSolver solver = new RiddersSolver();
        double min, max, expected, result, tolerance;

        min = 3.0; max = 4.0; expected = FastMath.PI;
        tolerance = FastMath.max(solver.getAbsoluteAccuracy(),
                    FastMath.abs(expected * solver.getRelativeAccuracy()));
        result = solver.solve(100, f, min, max);
        Assert.assertEquals(expected, result, tolerance);

        min = -1.0; max = 1.5; expected = 0.0;
        tolerance = FastMath.max(solver.getAbsoluteAccuracy(),
                    FastMath.abs(expected * solver.getRelativeAccuracy()));
        result = solver.solve(100, f, min, max);
        Assert.assertEquals(expected, result, tolerance);
    }

[202, 'src/test/java', 'org.apache.commons.math3.ml.distance', 'CanberraDistanceTest', 'testZero', 28, 32]

[202, 'src/test/java', 'org.apache.commons.math3.ml.distance', 'ChebyshevDistanceTest', 'testZero', 28, 32]

[202, 'src/test/java', 'org.apache.commons.math3.ml.distance', 'EuclideanDistanceTest', 'testZero', 29, 33]

[202, 'src/test/java', 'org.apache.commons.math3.ml.distance', 'ManhattanDistanceTest', 'testZero', 28, 32]

    @Test
    public void testZero() {
        final double[] a = { 0, 1, -2, 3.4, 5, -6.7, 89 };
        Assert.assertEquals(0, distance.compute(a, a), 0d);
    }

    @Test
    public void testZero() {
        final double[] a = { 0, 1, -2, 3.4, 5, -6.7, 89 };
        Assert.assertEquals(0, distance.compute(a, a), 0d);
    }

    @Test
    public void testZero() {
        final double[] a = { 0, 1, -2, 3.4, 5, -6.7, 89 };
        Assert.assertEquals(0, distance.compute(a, a), 0d);
    }

    @Test
    public void testZero() {
        final double[] a = { 0, 1, -2, 3.4, 5, -6.7, 89 };
        Assert.assertEquals(0, distance.compute(a, a), 0d);
    }

[218, 'src/test/java', 'org.apache.commons.math3.util', 'MathArraysTest', 'testEbeAddPrecondition', 75, 78]

[218, 'src/test/java', 'org.apache.commons.math3.util', 'MathArraysTest', 'testEbeSubtractPrecondition', 79, 82]

[218, 'src/test/java', 'org.apache.commons.math3.util', 'MathArraysTest', 'testEbeMultiplyPrecondition', 83, 86]

[218, 'src/test/java', 'org.apache.commons.math3.util', 'MathArraysTest', 'testEbeDividePrecondition', 87, 90]

    @Test(expected=DimensionMismatchException.class)
    public void testEbeAddPrecondition() {
        MathArrays.ebeAdd(new double[3], new double[4]);
    }

    @Test(expected=DimensionMismatchException.class)
    public void testEbeSubtractPrecondition() {
        MathArrays.ebeSubtract(new double[3], new double[4]);
    }

    @Test(expected=DimensionMismatchException.class)
    public void testEbeMultiplyPrecondition() {
        MathArrays.ebeMultiply(new double[3], new double[4]);
    }

    @Test(expected=DimensionMismatchException.class)
    public void testEbeDividePrecondition() {
        MathArrays.ebeDivide(new double[3], new double[4]);
    }

[220, 'src/test/java', 'org.apache.commons.math3.analysis.function', 'GaussianTest', 'testParametricUsage1', 95, 99]

[220, 'src/test/java', 'org.apache.commons.math3.analysis.function', 'GaussianTest', 'testParametricUsage4', 113, 117]

[220, 'src/test/java', 'org.apache.commons.math3.analysis.function', 'HarmonicOscillatorTest', 'testParametricUsage1', 85, 89]

[220, 'src/test/java', 'org.apache.commons.math3.analysis.function', 'HarmonicOscillatorTest', 'testParametricUsage3', 97, 101]

[220, 'src/test/java', 'org.apache.commons.math3.analysis.function', 'LogisticTest', 'testParametricUsage1', 106, 110]

[220, 'src/test/java', 'org.apache.commons.math3.analysis.function', 'LogisticTest', 'testParametricUsage3', 118, 122]

[220, 'src/test/java', 'org.apache.commons.math3.analysis.function', 'LogitTest', 'testParametricUsage1', 108, 112]

[220, 'src/test/java', 'org.apache.commons.math3.analysis.function', 'LogitTest', 'testParametricUsage3', 120, 124]

[220, 'src/test/java', 'org.apache.commons.math3.analysis.function', 'SigmoidTest', 'testParametricUsage1', 76, 80]

[220, 'src/test/java', 'org.apache.commons.math3.analysis.function', 'SigmoidTest', 'testParametricUsage3', 88, 92]

    @Test(expected=NullArgumentException.class)
    public void testParametricUsage1() {
        final Gaussian.Parametric g = new Gaussian.Parametric();
        g.value(0, null);
    }

    @Test(expected=NullArgumentException.class)
    public void testParametricUsage4() {
        final Gaussian.Parametric g = new Gaussian.Parametric();
        g.gradient(0, null);
    }

    @Test(expected=NullArgumentException.class)
    public void testParametricUsage1() {
        final HarmonicOscillator.Parametric g = new HarmonicOscillator.Parametric();
        g.value(0, null);
    }

    @Test(expected=NullArgumentException.class)
    public void testParametricUsage3() {
        final HarmonicOscillator.Parametric g = new HarmonicOscillator.Parametric();
        g.gradient(0, null);
    }

    @Test(expected=NullArgumentException.class)
    public void testParametricUsage1() {
        final Logistic.Parametric g = new Logistic.Parametric();
        g.value(0, null);
    }

    @Test(expected=NullArgumentException.class)
    public void testParametricUsage3() {
        final Logistic.Parametric g = new Logistic.Parametric();
        g.gradient(0, null);
    }

    @Test(expected=NullArgumentException.class)
    public void testParametricUsage1() {
        final Logit.Parametric g = new Logit.Parametric();
        g.value(0, null);
    }

    @Test(expected=NullArgumentException.class)
    public void testParametricUsage3() {
        final Logit.Parametric g = new Logit.Parametric();
        g.gradient(0, null);
    }

    @Test(expected=NullArgumentException.class)
    public void testParametricUsage1() {
        final Sigmoid.Parametric g = new Sigmoid.Parametric();
        g.value(0, null);
    }

    @Test(expected=NullArgumentException.class)
    public void testParametricUsage3() {
        final Sigmoid.Parametric g = new Sigmoid.Parametric();
        g.gradient(0, null);
    }

[233, 'src/test/java', 'org.apache.commons.math3.geometry.euclidean.threed', 'FieldRotationDSTest', 'testApplyToRotation', 781, 807]

[233, 'src/test/java', 'org.apache.commons.math3.geometry.euclidean.threed', 'FieldRotationDSTest', 'testApplyInverseToRotation', 867, 893]

[233, 'src/test/java', 'org.apache.commons.math3.geometry.euclidean.threed', 'FieldRotationDfpTest', 'testApplyToRotation', 618, 644]

[233, 'src/test/java', 'org.apache.commons.math3.geometry.euclidean.threed', 'FieldRotationDfpTest', 'testApplyInverseToRotation', 706, 732]

    @Test
    public void testApplyToRotation() throws MathIllegalArgumentException {

        FieldRotation<DerivativeStructure> r1       = new FieldRotation<DerivativeStructure>(createVector(2, -3, 5),
                                                                                             createAngle(1.7),
                                                                                             RotationConvention.VECTOR_OPERATOR);
        FieldRotation<DerivativeStructure> r2       = new FieldRotation<DerivativeStructure>(createVector(-1, 3, 2),
                                                                                             createAngle(0.3),
                                                                                             RotationConvention.VECTOR_OPERATOR);
        FieldRotation<DerivativeStructure> r3       = r2.applyTo(r1);
        FieldRotation<DerivativeStructure> r3Double = r2.applyTo(new Rotation(r1.getQ0().getReal(),
                                                                              r1.getQ1().getReal(),
                                                                              r1.getQ2().getReal(),
                                                                              r1.getQ3().getReal(),
                                                                              false));

        for (double x = -0.9; x < 0.9; x += 0.2) {
            for (double y = -0.9; y < 0.9; y += 0.2) {
                for (double z = -0.9; z < 0.9; z += 0.2) {
                    FieldVector3D<DerivativeStructure> u = createVector(x, y, z);
                    checkVector(r2.applyTo(r1.applyTo(u)), r3.applyTo(u));
                    checkVector(r2.applyTo(r1.applyTo(u)), r3Double.applyTo(u));
                }
            }
        }

    }

    @Test
    public void testApplyInverseToRotation() throws MathIllegalArgumentException {

        FieldRotation<DerivativeStructure> r1 = new FieldRotation<DerivativeStructure>(createVector(2, -3, 5),
                                                                                       createAngle(1.7),
                                                                                       RotationConvention.VECTOR_OPERATOR);
        FieldRotation<DerivativeStructure> r2 = new FieldRotation<DerivativeStructure>(createVector(-1, 3, 2),
                                                                                       createAngle(0.3),
                                                                                       RotationConvention.VECTOR_OPERATOR);
        FieldRotation<DerivativeStructure> r3 = r2.applyInverseTo(r1);
        FieldRotation<DerivativeStructure> r3Double = r2.applyInverseTo(new Rotation(r1.getQ0().getReal(),
                                                                                     r1.getQ1().getReal(),
                                                                                     r1.getQ2().getReal(),
                                                                                     r1.getQ3().getReal(),
                                                                                    false));

        for (double x = -0.9; x < 0.9; x += 0.2) {
            for (double y = -0.9; y < 0.9; y += 0.2) {
                for (double z = -0.9; z < 0.9; z += 0.2) {
                    FieldVector3D<DerivativeStructure> u = createVector(x, y, z);
                    checkVector(r2.applyInverseTo(r1.applyTo(u)), r3.applyTo(u));
                    checkVector(r2.applyInverseTo(r1.applyTo(u)), r3Double.applyTo(u));
                }
            }
        }

    }

    @Test
    public void testApplyToRotation() throws MathIllegalArgumentException {

        FieldRotation<Dfp> r1       = new FieldRotation<Dfp>(createVector(2, -3, 5),
                                                             createAngle(1.7),
                                                             RotationConvention.VECTOR_OPERATOR);
        FieldRotation<Dfp> r2       = new FieldRotation<Dfp>(createVector(-1, 3, 2),
                                                             createAngle(0.3),
                                                             RotationConvention.VECTOR_OPERATOR);
        FieldRotation<Dfp> r3       = r2.applyTo(r1);
        FieldRotation<Dfp> r3Double = r2.applyTo(new Rotation(r1.getQ0().getReal(),
                                                      r1.getQ1().getReal(),
                                                      r1.getQ2().getReal(),
                                                      r1.getQ3().getReal(),
                                                      false));

        for (double x = -0.9; x < 0.9; x += 0.2) {
            for (double y = -0.9; y < 0.9; y += 0.2) {
                for (double z = -0.9; z < 0.9; z += 0.2) {
                    FieldVector3D<Dfp> u = createVector(x, y, z);
                    checkVector(r2.applyTo(r1.applyTo(u)), r3.applyTo(u));
                    checkVector(r2.applyTo(r1.applyTo(u)), r3Double.applyTo(u));
                }
            }
        }

    }

    @Test
    public void testApplyInverseToRotation() throws MathIllegalArgumentException {

        FieldRotation<Dfp> r1 = new FieldRotation<Dfp>(createVector(2, -3, 5),
                                                       createAngle(1.7),
                                                       RotationConvention.VECTOR_OPERATOR);
        FieldRotation<Dfp> r2 = new FieldRotation<Dfp>(createVector(-1, 3, 2),
                                                       createAngle(0.3),
                                                       RotationConvention.VECTOR_OPERATOR);
        FieldRotation<Dfp> r3 = r2.applyInverseTo(r1);
        FieldRotation<Dfp> r3Double = r2.applyInverseTo(new Rotation(r1.getQ0().getReal(),
                                                             r1.getQ1().getReal(),
                                                             r1.getQ2().getReal(),
                                                             r1.getQ3().getReal(),
                                                             false));

        for (double x = -0.9; x < 0.9; x += 0.2) {
            for (double y = -0.9; y < 0.9; y += 0.2) {
                for (double z = -0.9; z < 0.9; z += 0.2) {
                    FieldVector3D<Dfp> u = createVector(x, y, z);
                    checkVector(r2.applyInverseTo(r1.applyTo(u)), r3.applyTo(u));
                    checkVector(r2.applyInverseTo(r1.applyTo(u)), r3Double.applyTo(u));
                }
            }
        }

    }

[250, 'src/test/java', 'org.apache.commons.math3.ode.nonstiff', 'DormandPrince54StepInterpolatorTest', 'checkClone', 110, 151]

[250, 'src/test/java', 'org.apache.commons.math3.ode.nonstiff', 'DormandPrince853StepInterpolatorTest', 'checklone', 110, 151]

[250, 'src/test/java', 'org.apache.commons.math3.ode.nonstiff', 'GraggBulirschStoerStepInterpolatorTest', 'checklone', 112, 153]

[250, 'src/test/java', 'org.apache.commons.math3.ode.nonstiff', 'HighamHall54StepInterpolatorTest', 'checkClone', 110, 151]

  @Test
  public void checkClone()
      throws DimensionMismatchException, NumberIsTooSmallException,
             MaxCountExceededException, NoBracketingException {
      TestProblem3 pb = new TestProblem3(0.9);
      double minStep = 0;
      double maxStep = pb.getFinalTime() - pb.getInitialTime();
      double scalAbsoluteTolerance = 1.0e-8;
      double scalRelativeTolerance = scalAbsoluteTolerance;
      DormandPrince54Integrator integ = new DormandPrince54Integrator(minStep, maxStep,
                                                                      scalAbsoluteTolerance,
                                                                      scalRelativeTolerance);
      integ.addStepHandler(new StepHandler() {
          public void handleStep(StepInterpolator interpolator, boolean isLast)
              throws MaxCountExceededException {
              StepInterpolator cloned = interpolator.copy();
              double tA = cloned.getPreviousTime();
              double tB = cloned.getCurrentTime();
              double halfStep = FastMath.abs(tB - tA) / 2;
              Assert.assertEquals(interpolator.getPreviousTime(), tA, 1.0e-12);
              Assert.assertEquals(interpolator.getCurrentTime(), tB, 1.0e-12);
              for (int i = 0; i < 10; ++i) {
                  double t = (i * tB + (9 - i) * tA) / 9;
                  interpolator.setInterpolatedTime(t);
                  Assert.assertTrue(FastMath.abs(cloned.getInterpolatedTime() - t) > (halfStep / 10));
                  cloned.setInterpolatedTime(t);
                  Assert.assertEquals(t, cloned.getInterpolatedTime(), 1.0e-12);
                  double[] referenceState = interpolator.getInterpolatedState();
                  double[] cloneState     = cloned.getInterpolatedState();
                  for (int j = 0; j < referenceState.length; ++j) {
                      Assert.assertEquals(referenceState[j], cloneState[j], 1.0e-12);
                  }
              }
          }
          public void init(double t0, double[] y0, double t) {
          }
      });
      integ.integrate(pb,
              pb.getInitialTime(), pb.getInitialState(),
              pb.getFinalTime(), new double[pb.getDimension()]);

  }

  @Test
  public void checklone()
      throws DimensionMismatchException, NumberIsTooSmallException,
             MaxCountExceededException, NoBracketingException {
    TestProblem3 pb = new TestProblem3(0.9);
    double minStep = 0;
    double maxStep = pb.getFinalTime() - pb.getInitialTime();
    double scalAbsoluteTolerance = 1.0e-8;
    double scalRelativeTolerance = scalAbsoluteTolerance;
    DormandPrince853Integrator integ = new DormandPrince853Integrator(minStep, maxStep,
                                                                      scalAbsoluteTolerance,
                                                                      scalRelativeTolerance);
    integ.addStepHandler(new StepHandler() {
        public void handleStep(StepInterpolator interpolator, boolean isLast)
            throws MaxCountExceededException {
            StepInterpolator cloned = interpolator.copy();
            double tA = cloned.getPreviousTime();
            double tB = cloned.getCurrentTime();
            double halfStep = FastMath.abs(tB - tA) / 2;
            Assert.assertEquals(interpolator.getPreviousTime(), tA, 1.0e-12);
            Assert.assertEquals(interpolator.getCurrentTime(), tB, 1.0e-12);
            for (int i = 0; i < 10; ++i) {
                double t = (i * tB + (9 - i) * tA) / 9;
                interpolator.setInterpolatedTime(t);
                Assert.assertTrue(FastMath.abs(cloned.getInterpolatedTime() - t) > (halfStep / 10));
                cloned.setInterpolatedTime(t);
                Assert.assertEquals(t, cloned.getInterpolatedTime(), 1.0e-12);
                double[] referenceState = interpolator.getInterpolatedState();
                double[] cloneState     = cloned.getInterpolatedState();
                for (int j = 0; j < referenceState.length; ++j) {
                    Assert.assertEquals(referenceState[j], cloneState[j], 1.0e-12);
                }
            }
        }
        public void init(double t0, double[] y0, double t) {
        }
    });
    integ.integrate(pb,
            pb.getInitialTime(), pb.getInitialState(),
            pb.getFinalTime(), new double[pb.getDimension()]);

  }

  @Test
  public void checklone()
      throws DimensionMismatchException, NumberIsTooSmallException,
             MaxCountExceededException, NoBracketingException {
    TestProblem3 pb = new TestProblem3(0.9);
    double minStep = 0;
    double maxStep = pb.getFinalTime() - pb.getInitialTime();
    double scalAbsoluteTolerance = 1.0e-8;
    double scalRelativeTolerance = scalAbsoluteTolerance;
    GraggBulirschStoerIntegrator integ = new GraggBulirschStoerIntegrator(minStep, maxStep,
                                                                          scalAbsoluteTolerance,
                                                                          scalRelativeTolerance);
    integ.addStepHandler(new StepHandler() {
        public void handleStep(StepInterpolator interpolator, boolean isLast)
            throws MaxCountExceededException {
            StepInterpolator cloned = interpolator.copy();
            double tA = cloned.getPreviousTime();
            double tB = cloned.getCurrentTime();
            double halfStep = FastMath.abs(tB - tA) / 2;
            Assert.assertEquals(interpolator.getPreviousTime(), tA, 1.0e-12);
            Assert.assertEquals(interpolator.getCurrentTime(), tB, 1.0e-12);
            for (int i = 0; i < 10; ++i) {
                double t = (i * tB + (9 - i) * tA) / 9;
                interpolator.setInterpolatedTime(t);
                Assert.assertTrue(FastMath.abs(cloned.getInterpolatedTime() - t) > (halfStep / 10));
                cloned.setInterpolatedTime(t);
                Assert.assertEquals(t, cloned.getInterpolatedTime(), 1.0e-12);
                double[] referenceState = interpolator.getInterpolatedState();
                double[] cloneState     = cloned.getInterpolatedState();
                for (int j = 0; j < referenceState.length; ++j) {
                    Assert.assertEquals(referenceState[j], cloneState[j], 1.0e-12);
                }
            }
        }
        public void init(double t0, double[] y0, double t) {
        }
    });
    integ.integrate(pb,
            pb.getInitialTime(), pb.getInitialState(),
            pb.getFinalTime(), new double[pb.getDimension()]);

  }

  @Test
  public void checkClone()
      throws DimensionMismatchException, NumberIsTooSmallException,
             MaxCountExceededException, NoBracketingException {
    TestProblem3 pb = new TestProblem3(0.9);
    double minStep = 0;
    double maxStep = pb.getFinalTime() - pb.getInitialTime();
    double scalAbsoluteTolerance = 1.0e-8;
    double scalRelativeTolerance = scalAbsoluteTolerance;
    HighamHall54Integrator integ = new HighamHall54Integrator(minStep, maxStep,
                                                              scalAbsoluteTolerance,
                                                              scalRelativeTolerance);
    integ.addStepHandler(new StepHandler() {
        public void handleStep(StepInterpolator interpolator, boolean isLast)
            throws MaxCountExceededException {
            StepInterpolator cloned = interpolator.copy();
            double tA = cloned.getPreviousTime();
            double tB = cloned.getCurrentTime();
            double halfStep = FastMath.abs(tB - tA) / 2;
            Assert.assertEquals(interpolator.getPreviousTime(), tA, 1.0e-12);
            Assert.assertEquals(interpolator.getCurrentTime(), tB, 1.0e-12);
            for (int i = 0; i < 10; ++i) {
                double t = (i * tB + (9 - i) * tA) / 9;
                interpolator.setInterpolatedTime(t);
                Assert.assertTrue(FastMath.abs(cloned.getInterpolatedTime() - t) > (halfStep / 10));
                cloned.setInterpolatedTime(t);
                Assert.assertEquals(t, cloned.getInterpolatedTime(), 1.0e-12);
                double[] referenceState = interpolator.getInterpolatedState();
                double[] cloneState     = cloned.getInterpolatedState();
                for (int j = 0; j < referenceState.length; ++j) {
                    Assert.assertEquals(referenceState[j], cloneState[j], 1.0e-12);
                }
            }
        }
        public void init(double t0, double[] y0, double t) {
        }
    });
    integ.integrate(pb,
            pb.getInitialTime(), pb.getInitialState(),
            pb.getFinalTime(), new double[pb.getDimension()]);

  }

[294, 'src/test/java', 'org.apache.commons.math3.analysis.integration', 'RombergIntegratorTest', 'testQuinticFunction', 66, 92]

[294, 'src/test/java', 'org.apache.commons.math3.analysis.integration', 'SimpsonIntegratorTest', 'testQuinticFunction', 65, 91]

[294, 'src/test/java', 'org.apache.commons.math3.analysis.integration', 'TrapezoidIntegratorTest', 'testQuinticFunction', 65, 92]

    @Test
    public void testQuinticFunction() {
        UnivariateFunction f = new QuinticFunction();
        UnivariateIntegrator integrator = new RombergIntegrator();
        double min, max, expected, result, tolerance;

        min = 0; max = 1; expected = -1.0/48;
        tolerance = FastMath.abs(expected * integrator.getRelativeAccuracy());
        result = integrator.integrate(100, f, min, max);
        Assert.assertTrue(integrator.getEvaluations() < 10);
        Assert.assertTrue(integrator.getIterations()  < 5);
        Assert.assertEquals(expected, result, tolerance);

        min = 0; max = 0.5; expected = 11.0/768;
        tolerance = FastMath.abs(expected * integrator.getRelativeAccuracy());
        result = integrator.integrate(100, f, min, max);
        Assert.assertTrue(integrator.getEvaluations() < 10);
        Assert.assertTrue(integrator.getIterations()  < 5);
        Assert.assertEquals(expected, result, tolerance);

        min = -1; max = 4; expected = 2048/3.0 - 78 + 1.0/48;
        tolerance = FastMath.abs(expected * integrator.getRelativeAccuracy());
        result = integrator.integrate(100, f, min, max);
        Assert.assertTrue(integrator.getEvaluations() < 10);
        Assert.assertTrue(integrator.getIterations()  < 5);
        Assert.assertEquals(expected, result, tolerance);
    }

    @Test
    public void testQuinticFunction() {
        UnivariateFunction f = new QuinticFunction();
        UnivariateIntegrator integrator = new SimpsonIntegrator();
        double min, max, expected, result, tolerance;

        min = 0; max = 1; expected = -1.0/48;
        tolerance = FastMath.abs(expected * integrator.getRelativeAccuracy());
        result = integrator.integrate(1000, f, min, max);
        Assert.assertTrue(integrator.getEvaluations() < 150);
        Assert.assertTrue(integrator.getIterations()  < 10);
        Assert.assertEquals(expected, result, tolerance);

        min = 0; max = 0.5; expected = 11.0/768;
        tolerance = FastMath.abs(expected * integrator.getRelativeAccuracy());
        result = integrator.integrate(1000, f, min, max);
        Assert.assertTrue(integrator.getEvaluations() < 100);
        Assert.assertTrue(integrator.getIterations()  < 10);
        Assert.assertEquals(expected, result, tolerance);

        min = -1; max = 4; expected = 2048/3.0 - 78 + 1.0/48;
        tolerance = FastMath.abs(expected * integrator.getRelativeAccuracy());
        result = integrator.integrate(1000, f, min, max);
        Assert.assertTrue(integrator.getEvaluations() < 150);
        Assert.assertTrue(integrator.getIterations()  < 10);
        Assert.assertEquals(expected, result, tolerance);
    }

    @Test
    public void testQuinticFunction() {
        UnivariateFunction f = new QuinticFunction();
        UnivariateIntegrator integrator = new TrapezoidIntegrator();
        double min, max, expected, result, tolerance;

        min = 0; max = 1; expected = -1.0/48;
        tolerance = FastMath.abs(expected * integrator.getRelativeAccuracy());
        result = integrator.integrate(10000, f, min, max);
        Assert.assertTrue(integrator.getEvaluations() < 5000);
        Assert.assertTrue(integrator.getIterations()  < 15);
        Assert.assertEquals(expected, result, tolerance);

        min = 0; max = 0.5; expected = 11.0/768;
        tolerance = FastMath.abs(expected * integrator.getRelativeAccuracy());
        result = integrator.integrate(10000, f, min, max);
        Assert.assertTrue(integrator.getEvaluations() < 2500);
        Assert.assertTrue(integrator.getIterations()  < 15);
        Assert.assertEquals(expected, result, tolerance);

        min = -1; max = 4; expected = 2048/3.0 - 78 + 1.0/48;
        tolerance = FastMath.abs(expected * integrator.getRelativeAccuracy());
        result = integrator.integrate(10000, f, min, max);
        Assert.assertTrue(integrator.getEvaluations() < 5000);
        Assert.assertTrue(integrator.getIterations()  < 15);
        Assert.assertEquals(expected, result, tolerance);

    }

[303, 'src/test/java', 'org.apache.commons.math3.ode.nonstiff', 'ClassicalRungeKuttaIntegratorTest', 'testBigStep', 204, 222]

[303, 'src/test/java', 'org.apache.commons.math3.ode.nonstiff', 'EulerIntegratorTest', 'testBigStep', 125, 144]

[303, 'src/test/java', 'org.apache.commons.math3.ode.nonstiff', 'GillIntegratorTest', 'testBigStep', 124, 142]

[303, 'src/test/java', 'org.apache.commons.math3.ode.nonstiff', 'LutherIntegratorTest', 'testBigStep', 204, 222]

[303, 'src/test/java', 'org.apache.commons.math3.ode.nonstiff', 'MidpointIntegratorTest', 'testBigStep', 125, 144]

[303, 'src/test/java', 'org.apache.commons.math3.ode.nonstiff', 'ThreeEighthesIntegratorTest', 'testBigStep', 124, 142]

  @Test
  public void testBigStep()
      throws DimensionMismatchException, NumberIsTooSmallException,
             MaxCountExceededException, NoBracketingException {

    TestProblem1 pb = new TestProblem1();
    double step = (pb.getFinalTime() - pb.getInitialTime()) * 0.2;

    FirstOrderIntegrator integ = new ClassicalRungeKuttaIntegrator(step);
    TestProblemHandler handler = new TestProblemHandler(pb, integ);
    integ.addStepHandler(handler);
    integ.integrate(pb, pb.getInitialTime(), pb.getInitialState(),
                    pb.getFinalTime(), new double[pb.getDimension()]);

    Assert.assertTrue(handler.getLastError() > 0.0004);
    Assert.assertTrue(handler.getMaximalValueError() > 0.005);
    Assert.assertEquals(0, handler.getMaximalTimeError(), 1.0e-12);

  }

  @Test
  public void testBigStep()
      throws DimensionMismatchException, NumberIsTooSmallException,
             MaxCountExceededException, NoBracketingException {

    TestProblem1 pb  = new TestProblem1();
    double step = (pb.getFinalTime() - pb.getInitialTime()) * 0.2;

    FirstOrderIntegrator integ = new EulerIntegrator(step);
    TestProblemHandler handler = new TestProblemHandler(pb, integ);
    integ.addStepHandler(handler);
    integ.integrate(pb,
                    pb.getInitialTime(), pb.getInitialState(),
                    pb.getFinalTime(), new double[pb.getDimension()]);

    Assert.assertTrue(handler.getLastError() > 0.01);
    Assert.assertTrue(handler.getMaximalValueError() > 0.2);
    Assert.assertEquals(0, handler.getMaximalTimeError(), 1.0e-12);

  }

  @Test
  public void testBigStep()
      throws DimensionMismatchException, NumberIsTooSmallException,
             MaxCountExceededException, NoBracketingException {

    TestProblem1 pb = new TestProblem1();
    double step = (pb.getFinalTime() - pb.getInitialTime()) * 0.2;

    FirstOrderIntegrator integ = new GillIntegrator(step);
    TestProblemHandler handler = new TestProblemHandler(pb, integ);
    integ.addStepHandler(handler);
    integ.integrate(pb, pb.getInitialTime(), pb.getInitialState(),
                    pb.getFinalTime(), new double[pb.getDimension()]);

    Assert.assertTrue(handler.getLastError() > 0.0004);
    Assert.assertTrue(handler.getMaximalValueError() > 0.005);
    Assert.assertEquals(0, handler.getMaximalTimeError(), 1.0e-12);

  }

    @Test
    public void testBigStep()
            throws DimensionMismatchException, NumberIsTooSmallException,
            MaxCountExceededException, NoBracketingException {

        TestProblem1 pb = new TestProblem1();
        double step = (pb.getFinalTime() - pb.getInitialTime()) * 0.2;

        FirstOrderIntegrator integ = new LutherIntegrator(step);
        TestProblemHandler handler = new TestProblemHandler(pb, integ);
        integ.addStepHandler(handler);
        integ.integrate(pb, pb.getInitialTime(), pb.getInitialState(),
                        pb.getFinalTime(), new double[pb.getDimension()]);

        Assert.assertTrue(handler.getLastError() > 0.00002);
        Assert.assertTrue(handler.getMaximalValueError() > 0.001);
        Assert.assertEquals(0, handler.getMaximalTimeError(), 1.0e-12);

    }

  @Test
  public void testBigStep()
      throws DimensionMismatchException, NumberIsTooSmallException,
             MaxCountExceededException, NoBracketingException {

    TestProblem1 pb  = new TestProblem1();
    double step = (pb.getFinalTime() - pb.getInitialTime()) * 0.2;

    FirstOrderIntegrator integ = new MidpointIntegrator(step);
    TestProblemHandler handler = new TestProblemHandler(pb, integ);
    integ.addStepHandler(handler);
    integ.integrate(pb,
                    pb.getInitialTime(), pb.getInitialState(),
                    pb.getFinalTime(), new double[pb.getDimension()]);

    Assert.assertTrue(handler.getLastError() > 0.01);
    Assert.assertTrue(handler.getMaximalValueError() > 0.05);
    Assert.assertEquals(0, handler.getMaximalTimeError(), 1.0e-12);

  }

  @Test
  public void testBigStep()
      throws DimensionMismatchException, NumberIsTooSmallException,
             MaxCountExceededException, NoBracketingException {

    TestProblem1 pb = new TestProblem1();
    double step = (pb.getFinalTime() - pb.getInitialTime()) * 0.2;

    FirstOrderIntegrator integ = new ThreeEighthesIntegrator(step);
    TestProblemHandler handler = new TestProblemHandler(pb, integ);
    integ.addStepHandler(handler);
    integ.integrate(pb, pb.getInitialTime(), pb.getInitialState(),
                    pb.getFinalTime(), new double[pb.getDimension()]);

    Assert.assertTrue(handler.getLastError() > 0.0004);
    Assert.assertTrue(handler.getMaximalValueError() > 0.005);
    Assert.assertEquals(0, handler.getMaximalTimeError(), 1.0e-12);

  }

[304, 'src/test/java', 'org.apache.commons.math3.linear', 'BlockFieldMatrixTest', 'testGetRowMatrix', 741, 760]

[304, 'src/test/java', 'org.apache.commons.math3.linear', 'FieldMatrixImplTest', 'testGetRowMatrix', 581, 602]

[304, 'src/test/java', 'org.apache.commons.math3.linear', 'FieldMatrixImplTest', 'testGetColumnMatrix', 625, 646]

    @Test
    public void testGetRowMatrix() {
        FieldMatrix<Fraction> m     = new BlockFieldMatrix<Fraction>(subTestData);
        FieldMatrix<Fraction> mRow0 = new BlockFieldMatrix<Fraction>(subRow0);
        FieldMatrix<Fraction> mRow3 = new BlockFieldMatrix<Fraction>(subRow3);
        Assert.assertEquals("Row0", mRow0, m.getRowMatrix(0));
        Assert.assertEquals("Row3", mRow3, m.getRowMatrix(3));
        try {
            m.getRowMatrix(-1);
            Assert.fail("Expecting OutOfRangeException");
        } catch (OutOfRangeException ex) {
            // expected
        }
        try {
            m.getRowMatrix(4);
            Assert.fail("Expecting OutOfRangeException");
        } catch (OutOfRangeException ex) {
            // expected
        }
    }

    @Test
    public void testGetRowMatrix() {
        FieldMatrix<Fraction> m = new Array2DRowFieldMatrix<Fraction>(subTestData);
        FieldMatrix<Fraction> mRow0 = new Array2DRowFieldMatrix<Fraction>(subRow0);
        FieldMatrix<Fraction> mRow3 = new Array2DRowFieldMatrix<Fraction>(subRow3);
        Assert.assertEquals("Row0", mRow0,
                m.getRowMatrix(0));
        Assert.assertEquals("Row3", mRow3,
                m.getRowMatrix(3));
        try {
            m.getRowMatrix(-1);
            Assert.fail("Expecting OutOfRangeException");
        } catch (OutOfRangeException ex) {
            // expected
        }
        try {
            m.getRowMatrix(4);
            Assert.fail("Expecting OutOfRangeException");
        } catch (OutOfRangeException ex) {
            // expected
        }
    }

    @Test
    public void testGetColumnMatrix() {
        FieldMatrix<Fraction> m = new Array2DRowFieldMatrix<Fraction>(subTestData);
        FieldMatrix<Fraction> mColumn1 = new Array2DRowFieldMatrix<Fraction>(subColumn1);
        FieldMatrix<Fraction> mColumn3 = new Array2DRowFieldMatrix<Fraction>(subColumn3);
        Assert.assertEquals("Column1", mColumn1,
                m.getColumnMatrix(1));
        Assert.assertEquals("Column3", mColumn3,
                m.getColumnMatrix(3));
        try {
            m.getColumnMatrix(-1);
            Assert.fail("Expecting OutOfRangeException");
        } catch (OutOfRangeException ex) {
            // expected
        }
        try {
            m.getColumnMatrix(4);
            Assert.fail("Expecting OutOfRangeException");
        } catch (OutOfRangeException ex) {
            // expected
        }
    }

[317, 'src/test/java', 'org.apache.commons.math3.geometry.euclidean.threed', 'FieldRotationDSTest', 'testComposeVectorOperator', 809, 836]

[317, 'src/test/java', 'org.apache.commons.math3.geometry.euclidean.threed', 'FieldRotationDSTest', 'testComposeFrameTransform', 838, 865]

[317, 'src/test/java', 'org.apache.commons.math3.geometry.euclidean.threed', 'FieldRotationDSTest', 'testComposeInverseVectorOperator', 895, 922]

[317, 'src/test/java', 'org.apache.commons.math3.geometry.euclidean.threed', 'FieldRotationDSTest', 'testComposeInverseframeTransform', 924, 951]

[317, 'src/test/java', 'org.apache.commons.math3.geometry.euclidean.threed', 'FieldRotationDfpTest', 'testComposeVectorOperator', 646, 673]

[317, 'src/test/java', 'org.apache.commons.math3.geometry.euclidean.threed', 'FieldRotationDfpTest', 'testComposeInverseVectorOperator', 734, 761]

    @Test
    public void testComposeVectorOperator() throws MathIllegalArgumentException {

        FieldRotation<DerivativeStructure> r1       = new FieldRotation<DerivativeStructure>(createVector(2, -3, 5),
                                                                                             createAngle(1.7),
                                                                                             RotationConvention.VECTOR_OPERATOR);
        FieldRotation<DerivativeStructure> r2       = new FieldRotation<DerivativeStructure>(createVector(-1, 3, 2),
                                                                                             createAngle(0.3),
                                                                                             RotationConvention.VECTOR_OPERATOR);
        FieldRotation<DerivativeStructure> r3       = r2.compose(r1, RotationConvention.VECTOR_OPERATOR);
        FieldRotation<DerivativeStructure> r3Double = r2.compose(new Rotation(r1.getQ0().getReal(),
                                                                              r1.getQ1().getReal(),
                                                                              r1.getQ2().getReal(),
                                                                              r1.getQ3().getReal(),
                                                                              false),
                                                                 RotationConvention.VECTOR_OPERATOR);

        for (double x = -0.9; x < 0.9; x += 0.2) {
            for (double y = -0.9; y < 0.9; y += 0.2) {
                for (double z = -0.9; z < 0.9; z += 0.2) {
                    FieldVector3D<DerivativeStructure> u = createVector(x, y, z);
                    checkVector(r2.applyTo(r1.applyTo(u)), r3.applyTo(u));
                    checkVector(r2.applyTo(r1.applyTo(u)), r3Double.applyTo(u));
                }
            }
        }

    }

    @Test
    public void testComposeFrameTransform() throws MathIllegalArgumentException {

        FieldRotation<DerivativeStructure> r1       = new FieldRotation<DerivativeStructure>(createVector(2, -3, 5),
                                                                                             createAngle(1.7),
                                                                                             RotationConvention.FRAME_TRANSFORM);
        FieldRotation<DerivativeStructure> r2       = new FieldRotation<DerivativeStructure>(createVector(-1, 3, 2),
                                                                                             createAngle(0.3),
                                                                                             RotationConvention.FRAME_TRANSFORM);
        FieldRotation<DerivativeStructure> r3       = r2.compose(r1, RotationConvention.FRAME_TRANSFORM);
        FieldRotation<DerivativeStructure> r3Double = r2.compose(new Rotation(r1.getQ0().getReal(),
                                                                              r1.getQ1().getReal(),
                                                                              r1.getQ2().getReal(),
                                                                              r1.getQ3().getReal(),
                                                                              false),
                                                                 RotationConvention.FRAME_TRANSFORM);

        for (double x = -0.9; x < 0.9; x += 0.2) {
            for (double y = -0.9; y < 0.9; y += 0.2) {
                for (double z = -0.9; z < 0.9; z += 0.2) {
                    FieldVector3D<DerivativeStructure> u = createVector(x, y, z);
                    checkVector(r1.applyTo(r2.applyTo(u)), r3.applyTo(u));
                    checkVector(r1.applyTo(r2.applyTo(u)), r3Double.applyTo(u));
                }
            }
        }

    }

    @Test
    public void testComposeInverseVectorOperator() throws MathIllegalArgumentException {

        FieldRotation<DerivativeStructure> r1 = new FieldRotation<DerivativeStructure>(createVector(2, -3, 5),
                                                                                       createAngle(1.7),
                                                                                       RotationConvention.VECTOR_OPERATOR);
        FieldRotation<DerivativeStructure> r2 = new FieldRotation<DerivativeStructure>(createVector(-1, 3, 2),
                                                                                       createAngle(0.3),
                                                                                       RotationConvention.VECTOR_OPERATOR);
        FieldRotation<DerivativeStructure> r3 = r2.composeInverse(r1, RotationConvention.VECTOR_OPERATOR);
        FieldRotation<DerivativeStructure> r3Double = r2.composeInverse(new Rotation(r1.getQ0().getReal(),
                                                                                     r1.getQ1().getReal(),
                                                                                     r1.getQ2().getReal(),
                                                                                     r1.getQ3().getReal(),
                                                                                     false),
                                                                        RotationConvention.VECTOR_OPERATOR);

        for (double x = -0.9; x < 0.9; x += 0.2) {
            for (double y = -0.9; y < 0.9; y += 0.2) {
                for (double z = -0.9; z < 0.9; z += 0.2) {
                    FieldVector3D<DerivativeStructure> u = createVector(x, y, z);
                    checkVector(r2.applyInverseTo(r1.applyTo(u)), r3.applyTo(u));
                    checkVector(r2.applyInverseTo(r1.applyTo(u)), r3Double.applyTo(u));
                }
            }
        }

    }

    @Test
    public void testComposeInverseframeTransform() throws MathIllegalArgumentException {

        FieldRotation<DerivativeStructure> r1 = new FieldRotation<DerivativeStructure>(createVector(2, -3, 5),
                                                                                       createAngle(1.7),
                                                                                       RotationConvention.FRAME_TRANSFORM);
        FieldRotation<DerivativeStructure> r2 = new FieldRotation<DerivativeStructure>(createVector(-1, 3, 2),
                                                                                       createAngle(0.3),
                                                                                       RotationConvention.FRAME_TRANSFORM);
        FieldRotation<DerivativeStructure> r3 = r2.composeInverse(r1, RotationConvention.FRAME_TRANSFORM);
        FieldRotation<DerivativeStructure> r3Double = r2.composeInverse(new Rotation(r1.getQ0().getReal(),
                                                                                     r1.getQ1().getReal(),
                                                                                     r1.getQ2().getReal(),
                                                                                     r1.getQ3().getReal(),
                                                                                     false),
                                                                        RotationConvention.FRAME_TRANSFORM);

        for (double x = -0.9; x < 0.9; x += 0.2) {
            for (double y = -0.9; y < 0.9; y += 0.2) {
                for (double z = -0.9; z < 0.9; z += 0.2) {
                    FieldVector3D<DerivativeStructure> u = createVector(x, y, z);
                    checkVector(r1.applyTo(r2.applyInverseTo(u)), r3.applyTo(u));
                    checkVector(r1.applyTo(r2.applyInverseTo(u)), r3Double.applyTo(u));
                }
            }
        }

    }

    @Test
    public void testComposeVectorOperator() throws MathIllegalArgumentException {

        FieldRotation<Dfp> r1       = new FieldRotation<Dfp>(createVector(2, -3, 5),
                                                             createAngle(1.7),
                                                             RotationConvention.VECTOR_OPERATOR);
        FieldRotation<Dfp> r2       = new FieldRotation<Dfp>(createVector(-1, 3, 2),
                                                             createAngle(0.3),
                                                             RotationConvention.VECTOR_OPERATOR);
        FieldRotation<Dfp> r3       = r2.compose(r1, RotationConvention.VECTOR_OPERATOR);
        FieldRotation<Dfp> r3Double = r2.compose(new Rotation(r1.getQ0().getReal(),
                                                      r1.getQ1().getReal(),
                                                      r1.getQ2().getReal(),
                                                      r1.getQ3().getReal(),
                                                      false),
                                                 RotationConvention.VECTOR_OPERATOR);

        for (double x = -0.9; x < 0.9; x += 0.2) {
            for (double y = -0.9; y < 0.9; y += 0.2) {
                for (double z = -0.9; z < 0.9; z += 0.2) {
                    FieldVector3D<Dfp> u = createVector(x, y, z);
                    checkVector(r2.applyTo(r1.applyTo(u)), r3.applyTo(u));
                    checkVector(r2.applyTo(r1.applyTo(u)), r3Double.applyTo(u));
                }
            }
        }

    }

    @Test
    public void testComposeInverseVectorOperator() throws MathIllegalArgumentException {

        FieldRotation<Dfp> r1 = new FieldRotation<Dfp>(createVector(2, -3, 5),
                                                       createAngle(1.7),
                                                       RotationConvention.VECTOR_OPERATOR);
        FieldRotation<Dfp> r2 = new FieldRotation<Dfp>(createVector(-1, 3, 2),
                                                       createAngle(0.3),
                                                       RotationConvention.VECTOR_OPERATOR);
        FieldRotation<Dfp> r3 = r2.composeInverse(r1, RotationConvention.VECTOR_OPERATOR);
        FieldRotation<Dfp> r3Double = r2.composeInverse(new Rotation(r1.getQ0().getReal(),
                                                             r1.getQ1().getReal(),
                                                             r1.getQ2().getReal(),
                                                             r1.getQ3().getReal(),
                                                             false),
                                                        RotationConvention.VECTOR_OPERATOR);

        for (double x = -0.9; x < 0.9; x += 0.2) {
            for (double y = -0.9; y < 0.9; y += 0.2) {
                for (double z = -0.9; z < 0.9; z += 0.2) {
                    FieldVector3D<Dfp> u = createVector(x, y, z);
                    checkVector(r2.applyInverseTo(r1.applyTo(u)), r3.applyTo(u));
                    checkVector(r2.applyInverseTo(r1.applyTo(u)), r3Double.applyTo(u));
                }
            }
        }

    }

[331, 'src/test/java', 'org.apache.commons.math3.ode.nonstiff', 'DormandPrince54IntegratorTest', 'testBackward', 103, 126]

[331, 'src/test/java', 'org.apache.commons.math3.ode.nonstiff', 'DormandPrince853IntegratorTest', 'testBackward', 218, 241]

[331, 'src/test/java', 'org.apache.commons.math3.ode.nonstiff', 'GraggBulirschStoerIntegratorTest', 'testBackward', 88, 111]

[331, 'src/test/java', 'org.apache.commons.math3.ode.nonstiff', 'HighamHall54IntegratorTest', 'testBackward', 136, 159]

  @Test
  public void testBackward()
      throws DimensionMismatchException, NumberIsTooSmallException,
             MaxCountExceededException, NoBracketingException {

      TestProblem5 pb = new TestProblem5();
      double minStep = 0;
      double maxStep = pb.getFinalTime() - pb.getInitialTime();
      double scalAbsoluteTolerance = 1.0e-8;
      double scalRelativeTolerance = 0.01 * scalAbsoluteTolerance;

      FirstOrderIntegrator integ = new DormandPrince54Integrator(minStep, maxStep,
                                                                 scalAbsoluteTolerance,
                                                                 scalRelativeTolerance);
      TestProblemHandler handler = new TestProblemHandler(pb, integ);
      integ.addStepHandler(handler);
      integ.integrate(pb, pb.getInitialTime(), pb.getInitialState(),
                      pb.getFinalTime(), new double[pb.getDimension()]);

      Assert.assertTrue(handler.getLastError() < 2.0e-7);
      Assert.assertTrue(handler.getMaximalValueError() < 2.0e-7);
      Assert.assertEquals(0, handler.getMaximalTimeError(), 1.0e-12);
      Assert.assertEquals("Dormand-Prince 5(4)", integ.getName());
  }

  @Test
  public void testBackward()
      throws DimensionMismatchException, NumberIsTooSmallException,
             MaxCountExceededException, NoBracketingException {

      TestProblem5 pb = new TestProblem5();
      double minStep = 0;
      double maxStep = pb.getFinalTime() - pb.getInitialTime();
      double scalAbsoluteTolerance = 1.0e-8;
      double scalRelativeTolerance = 0.01 * scalAbsoluteTolerance;

      FirstOrderIntegrator integ = new DormandPrince853Integrator(minStep, maxStep,
                                                                  scalAbsoluteTolerance,
                                                                  scalRelativeTolerance);
      TestProblemHandler handler = new TestProblemHandler(pb, integ);
      integ.addStepHandler(handler);
      integ.integrate(pb, pb.getInitialTime(), pb.getInitialState(),
                      pb.getFinalTime(), new double[pb.getDimension()]);

      Assert.assertTrue(handler.getLastError() < 1.1e-7);
      Assert.assertTrue(handler.getMaximalValueError() < 1.1e-7);
      Assert.assertEquals(0, handler.getMaximalTimeError(), 1.0e-12);
      Assert.assertEquals("Dormand-Prince 8 (5, 3)", integ.getName());
  }

  @Test
  public void testBackward()
      throws DimensionMismatchException, NumberIsTooSmallException,
             MaxCountExceededException, NoBracketingException {

      TestProblem5 pb = new TestProblem5();
      double minStep = 0;
      double maxStep = pb.getFinalTime() - pb.getInitialTime();
      double scalAbsoluteTolerance = 1.0e-8;
      double scalRelativeTolerance = 0.01 * scalAbsoluteTolerance;

      FirstOrderIntegrator integ = new GraggBulirschStoerIntegrator(minStep, maxStep,
                                                                    scalAbsoluteTolerance,
                                                                    scalRelativeTolerance);
      TestProblemHandler handler = new TestProblemHandler(pb, integ);
      integ.addStepHandler(handler);
      integ.integrate(pb, pb.getInitialTime(), pb.getInitialState(),
                      pb.getFinalTime(), new double[pb.getDimension()]);

      Assert.assertTrue(handler.getLastError() < 7.5e-9);
      Assert.assertTrue(handler.getMaximalValueError() < 8.1e-9);
      Assert.assertEquals(0, handler.getMaximalTimeError(), 1.0e-12);
      Assert.assertEquals("Gragg-Bulirsch-Stoer", integ.getName());
  }

  @Test
  public void testBackward()
      throws DimensionMismatchException, NumberIsTooSmallException,
             MaxCountExceededException, NoBracketingException {

      TestProblem5 pb = new TestProblem5();
      double minStep = 0;
      double maxStep = pb.getFinalTime() - pb.getInitialTime();
      double scalAbsoluteTolerance = 1.0e-8;
      double scalRelativeTolerance = 0.01 * scalAbsoluteTolerance;

      FirstOrderIntegrator integ = new HighamHall54Integrator(minStep, maxStep,
                                                              scalAbsoluteTolerance,
                                                              scalRelativeTolerance);
      TestProblemHandler handler = new TestProblemHandler(pb, integ);
      integ.addStepHandler(handler);
      integ.integrate(pb, pb.getInitialTime(), pb.getInitialState(),
                      pb.getFinalTime(), new double[pb.getDimension()]);

      Assert.assertTrue(handler.getLastError() < 5.0e-7);
      Assert.assertTrue(handler.getMaximalValueError() < 5.0e-7);
      Assert.assertEquals(0, handler.getMaximalTimeError(), 1.0e-12);
      Assert.assertEquals("Higham-Hall 5(4)", integ.getName());
  }

[381, 'src/test/java', 'org.apache.commons.math3.optim.nonlinear.scalar.noderiv', 'SimplexOptimizerMultiDirectionalTest', 'testMinimize1', 48, 67]

[381, 'src/test/java', 'org.apache.commons.math3.optim.nonlinear.scalar.noderiv', 'SimplexOptimizerNelderMeadTest', 'testMinimize1', 53, 72]

[381, 'src/test/java', 'org.apache.commons.math3.optim.nonlinear.scalar.noderiv', 'SimplexOptimizerNelderMeadTest', 'testMaximize1', 95, 114]

    @Test
    public void testMinimize1() {
        SimplexOptimizer optimizer = new SimplexOptimizer(1e-11, 1e-30);
        final FourExtrema fourExtrema = new FourExtrema();

        final PointValuePair optimum
            = optimizer.optimize(new MaxEval(200),
                                 new ObjectiveFunction(fourExtrema),
                                 GoalType.MINIMIZE,
                                 new InitialGuess(new double[] { -3, 0 }),
                                 new MultiDirectionalSimplex(new double[] { 0.2, 0.2 }));
        Assert.assertEquals(fourExtrema.xM, optimum.getPoint()[0], 4e-6);
        Assert.assertEquals(fourExtrema.yP, optimum.getPoint()[1], 3e-6);
        Assert.assertEquals(fourExtrema.valueXmYp, optimum.getValue(), 8e-13);
        Assert.assertTrue(optimizer.getEvaluations() > 120);
        Assert.assertTrue(optimizer.getEvaluations() < 150);

        // Check that the number of iterations is updated (MATH-949).
        Assert.assertTrue(optimizer.getIterations() > 0);
    }

    @Test
    public void testMinimize1() {
        SimplexOptimizer optimizer = new SimplexOptimizer(1e-10, 1e-30);
        final FourExtrema fourExtrema = new FourExtrema();

        final PointValuePair optimum
            = optimizer.optimize(new MaxEval(100),
                                 new ObjectiveFunction(fourExtrema),
                                 GoalType.MINIMIZE,
                                 new InitialGuess(new double[] { -3, 0 }),
                                 new NelderMeadSimplex(new double[] { 0.2, 0.2 }));
        Assert.assertEquals(fourExtrema.xM, optimum.getPoint()[0], 2e-7);
        Assert.assertEquals(fourExtrema.yP, optimum.getPoint()[1], 2e-5);
        Assert.assertEquals(fourExtrema.valueXmYp, optimum.getValue(), 6e-12);
        Assert.assertTrue(optimizer.getEvaluations() > 60);
        Assert.assertTrue(optimizer.getEvaluations() < 90);

        // Check that the number of iterations is updated (MATH-949).
        Assert.assertTrue(optimizer.getIterations() > 0);
    }

    @Test
    public void testMaximize1() {
        SimplexOptimizer optimizer = new SimplexOptimizer(1e-10, 1e-30);
        final FourExtrema fourExtrema = new FourExtrema();

        final PointValuePair optimum
            = optimizer.optimize(new MaxEval(100),
                                 new ObjectiveFunction(fourExtrema),
                                 GoalType.MAXIMIZE,
                                 new InitialGuess(new double[] { -3, 0 }),
                                 new NelderMeadSimplex(new double[] { 0.2, 0.2 }));
        Assert.assertEquals(fourExtrema.xM, optimum.getPoint()[0], 1e-5);
        Assert.assertEquals(fourExtrema.yM, optimum.getPoint()[1], 3e-6);
        Assert.assertEquals(fourExtrema.valueXmYm, optimum.getValue(), 3e-12);
        Assert.assertTrue(optimizer.getEvaluations() > 60);
        Assert.assertTrue(optimizer.getEvaluations() < 90);

        // Check that the number of iterations is updated (MATH-949).
        Assert.assertTrue(optimizer.getIterations() > 0);
    }

[415, 'src/test/java', 'org.apache.commons.math3.geometry.euclidean.threed', 'FieldRotationDSTest', 'testAxisAngleVectorOperator', 88, 127]

[415, 'src/test/java', 'org.apache.commons.math3.geometry.euclidean.threed', 'FieldRotationDSTest', 'testAxisAngleFrameTransform', 129, 168]

[415, 'src/test/java', 'org.apache.commons.math3.geometry.euclidean.threed', 'FieldRotationDfpTest', 'testAxisAngleVectorOperator', 87, 126]

[415, 'src/test/java', 'org.apache.commons.math3.geometry.euclidean.threed', 'FieldRotationDfpTest', 'testAxisAngleFrameTransform', 128, 167]

    @Test
    public void testAxisAngleVectorOperator() throws MathIllegalArgumentException {

        FieldRotation<DerivativeStructure> r = new FieldRotation<DerivativeStructure>(createAxis(10, 10, 10),
                                                                                      createAngle(2 * FastMath.PI / 3) ,
                                                                                      RotationConvention.VECTOR_OPERATOR);
        checkVector(r.applyTo(createVector(1, 0, 0)), createVector(0, 1, 0));
        checkVector(r.applyTo(createVector(0, 1, 0)), createVector(0, 0, 1));
        checkVector(r.applyTo(createVector(0, 0, 1)), createVector(1, 0, 0));
        double s = 1 / FastMath.sqrt(3);
        checkVector(r.getAxis(RotationConvention.VECTOR_OPERATOR), createVector( s,  s,  s));
        checkVector(r.getAxis(RotationConvention.FRAME_TRANSFORM), createVector(-s, -s, -s));
        checkAngle(r.getAngle(), 2 * FastMath.PI / 3);

        try {
            new FieldRotation<DerivativeStructure>(createAxis(0, 0, 0),
                                                   createAngle(2 * FastMath.PI / 3),
                                                   RotationConvention.VECTOR_OPERATOR);
            Assert.fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException e) {
        }

        r = new FieldRotation<DerivativeStructure>(createAxis(0, 0, 1),
                                                   createAngle(1.5 * FastMath.PI),
                                                   RotationConvention.VECTOR_OPERATOR);
        checkVector(r.getAxis(RotationConvention.VECTOR_OPERATOR), createVector(0, 0, -1));
        checkVector(r.getAxis(RotationConvention.FRAME_TRANSFORM), createVector(0, 0, +1));
        checkAngle(r.getAngle(), 0.5 * FastMath.PI);

        r = new FieldRotation<DerivativeStructure>(createAxis(0, 1, 0),
                                                   createAngle(FastMath.PI),
                                                   RotationConvention.VECTOR_OPERATOR);
        checkVector(r.getAxis(RotationConvention.VECTOR_OPERATOR), createVector(0, +1, 0));
        checkVector(r.getAxis(RotationConvention.FRAME_TRANSFORM), createVector(0, -1, 0));
        checkAngle(r.getAngle(), FastMath.PI);

        checkVector(createRotation(1, 0, 0, 0, false).getAxis(RotationConvention.VECTOR_OPERATOR), createVector(+1, 0, 0));
        checkVector(createRotation(1, 0, 0, 0, false).getAxis(RotationConvention.FRAME_TRANSFORM), createVector(-1, 0, 0));

    }

    @Test
    public void testAxisAngleFrameTransform() throws MathIllegalArgumentException {

        FieldRotation<DerivativeStructure> r = new FieldRotation<DerivativeStructure>(createAxis(10, 10, 10),
                                                                                      createAngle(2 * FastMath.PI / 3) ,
                                                                                      RotationConvention.FRAME_TRANSFORM);
        checkVector(r.applyTo(createVector(1, 0, 0)), createVector(0, 0, 1));
        checkVector(r.applyTo(createVector(0, 1, 0)), createVector(1, 0, 0));
        checkVector(r.applyTo(createVector(0, 0, 1)), createVector(0, 1, 0));
        double s = 1 / FastMath.sqrt(3);
        checkVector(r.getAxis(RotationConvention.FRAME_TRANSFORM), createVector( s,  s,  s));
        checkVector(r.getAxis(RotationConvention.VECTOR_OPERATOR), createVector(-s, -s, -s));
        checkAngle(r.getAngle(), 2 * FastMath.PI / 3);

        try {
            new FieldRotation<DerivativeStructure>(createAxis(0, 0, 0),
                                                   createAngle(2 * FastMath.PI / 3),
                                                   RotationConvention.FRAME_TRANSFORM);
            Assert.fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException e) {
        }

        r = new FieldRotation<DerivativeStructure>(createAxis(0, 0, 1),
                                                   createAngle(1.5 * FastMath.PI),
                                                   RotationConvention.FRAME_TRANSFORM);
        checkVector(r.getAxis(RotationConvention.FRAME_TRANSFORM), createVector(0, 0, -1));
        checkVector(r.getAxis(RotationConvention.VECTOR_OPERATOR), createVector(0, 0, +1));
        checkAngle(r.getAngle(), 0.5 * FastMath.PI);

        r = new FieldRotation<DerivativeStructure>(createAxis(0, 1, 0),
                                                   createAngle(FastMath.PI),
                                                   RotationConvention.FRAME_TRANSFORM);
        checkVector(r.getAxis(RotationConvention.FRAME_TRANSFORM), createVector(0, +1, 0));
        checkVector(r.getAxis(RotationConvention.VECTOR_OPERATOR), createVector(0, -1, 0));
        checkAngle(r.getAngle(), FastMath.PI);

        checkVector(createRotation(1, 0, 0, 0, false).getAxis(RotationConvention.FRAME_TRANSFORM), createVector(-1, 0, 0));
        checkVector(createRotation(1, 0, 0, 0, false).getAxis(RotationConvention.VECTOR_OPERATOR), createVector(+1, 0, 0));

    }

    @Test
    public void testAxisAngleVectorOperator() throws MathIllegalArgumentException {

        FieldRotation<Dfp> r = new FieldRotation<Dfp>(createAxis(10, 10, 10),
                                                      createAngle(2 * FastMath.PI / 3) ,
                                                      RotationConvention.VECTOR_OPERATOR);
        checkVector(r.applyTo(createVector(1, 0, 0)), createVector(0, 1, 0));
        checkVector(r.applyTo(createVector(0, 1, 0)), createVector(0, 0, 1));
        checkVector(r.applyTo(createVector(0, 0, 1)), createVector(1, 0, 0));
        double s = 1 / FastMath.sqrt(3);
        checkVector(r.getAxis(RotationConvention.VECTOR_OPERATOR), createVector( s,  s,  s));
        checkVector(r.getAxis(RotationConvention.FRAME_TRANSFORM), createVector(-s, -s, -s));
        checkAngle(r.getAngle(), 2 * FastMath.PI / 3);

        try {
            new FieldRotation<Dfp>(createAxis(0, 0, 0),
                                   createAngle(2 * FastMath.PI / 3),
                                   RotationConvention.VECTOR_OPERATOR);
            Assert.fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException e) {
        }

        r = new FieldRotation<Dfp>(createAxis(0, 0, 1),
                                   createAngle(1.5 * FastMath.PI),
                                   RotationConvention.VECTOR_OPERATOR);
        checkVector(r.getAxis(RotationConvention.VECTOR_OPERATOR), createVector(0, 0, -1));
        checkVector(r.getAxis(RotationConvention.FRAME_TRANSFORM), createVector(0, 0, +1));
        checkAngle(r.getAngle(), 0.5 * FastMath.PI);

        r = new FieldRotation<Dfp>(createAxis(0, 1, 0),
                                   createAngle(FastMath.PI),
                                   RotationConvention.VECTOR_OPERATOR);
        checkVector(r.getAxis(RotationConvention.VECTOR_OPERATOR), createVector(0, +1, 0));
        checkVector(r.getAxis(RotationConvention.FRAME_TRANSFORM), createVector(0, -1, 0));
        checkAngle(r.getAngle(), FastMath.PI);

        checkVector(createRotation(1, 0, 0, 0, false).getAxis(RotationConvention.VECTOR_OPERATOR), createVector(+1, 0, 0));
        checkVector(createRotation(1, 0, 0, 0, false).getAxis(RotationConvention.FRAME_TRANSFORM), createVector(-1, 0, 0));

    }

    @Test
    public void testAxisAngleFrameTransform() throws MathIllegalArgumentException {

        FieldRotation<Dfp> r = new FieldRotation<Dfp>(createAxis(10, 10, 10),
                                                      createAngle(2 * FastMath.PI / 3) ,
                                                      RotationConvention.FRAME_TRANSFORM);
        checkVector(r.applyTo(createVector(1, 0, 0)), createVector(0, 0, 1));
        checkVector(r.applyTo(createVector(0, 1, 0)), createVector(1, 0, 0));
        checkVector(r.applyTo(createVector(0, 0, 1)), createVector(0, 1, 0));
        double s = 1 / FastMath.sqrt(3);
        checkVector(r.getAxis(RotationConvention.FRAME_TRANSFORM), createVector( s,  s,  s));
        checkVector(r.getAxis(RotationConvention.VECTOR_OPERATOR), createVector(-s, -s, -s));
        checkAngle(r.getAngle(), 2 * FastMath.PI / 3);

        try {
            new FieldRotation<Dfp>(createAxis(0, 0, 0),
                                   createAngle(2 * FastMath.PI / 3),
                                   RotationConvention.FRAME_TRANSFORM);
            Assert.fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException e) {
        }

        r = new FieldRotation<Dfp>(createAxis(0, 0, 1),
                                   createAngle(1.5 * FastMath.PI),
                                   RotationConvention.FRAME_TRANSFORM);
        checkVector(r.getAxis(RotationConvention.FRAME_TRANSFORM), createVector(0, 0, -1));
        checkVector(r.getAxis(RotationConvention.VECTOR_OPERATOR), createVector(0, 0, +1));
        checkAngle(r.getAngle(), 0.5 * FastMath.PI);

        r = new FieldRotation<Dfp>(createAxis(0, 1, 0),
                                   createAngle(FastMath.PI),
                                   RotationConvention.FRAME_TRANSFORM);
        checkVector(r.getAxis(RotationConvention.FRAME_TRANSFORM), createVector(0, +1, 0));
        checkVector(r.getAxis(RotationConvention.VECTOR_OPERATOR), createVector(0, -1, 0));
        checkAngle(r.getAngle(), FastMath.PI);

        checkVector(createRotation(1, 0, 0, 0, false).getAxis(RotationConvention.FRAME_TRANSFORM), createVector(-1, 0, 0));
        checkVector(createRotation(1, 0, 0, 0, false).getAxis(RotationConvention.VECTOR_OPERATOR), createVector(+1, 0, 0));

    }

[430, 'src/test/java', 'org.apache.commons.math3.random', 'RandomDataGeneratorTest', 'testNextIntIAE', 98, 106]

[430, 'src/test/java', 'org.apache.commons.math3.random', 'RandomDataGeneratorTest', 'testNextLongIAE', 173, 181]

[430, 'src/test/java', 'org.apache.commons.math3.random', 'RandomDataGeneratorTest', 'testNextSecureLongIAE', 248, 256]

[430, 'src/test/java', 'org.apache.commons.math3.random', 'RandomDataGeneratorTest', 'testNextSecureIntIAE', 305, 313]

    @Test
    public void testNextIntIAE() {
        try {
            randomData.nextInt(4, 3);
            Assert.fail("MathIllegalArgumentException expected");
        } catch (MathIllegalArgumentException ex) {
            // ignored
        }
    }

    @Test
    public void testNextLongIAE() {
        try {
            randomData.nextLong(4, 3);
            Assert.fail("MathIllegalArgumentException expected");
        } catch (MathIllegalArgumentException ex) {
            // ignored
        }
    }

    @Test
    public void testNextSecureLongIAE() {
        try {
            randomData.nextSecureLong(4, 3);
            Assert.fail("MathIllegalArgumentException expected");
        } catch (MathIllegalArgumentException ex) {
            // ignored
        }
    }

    @Test
    public void testNextSecureIntIAE() {
        try {
            randomData.nextSecureInt(4, 3);
            Assert.fail("MathIllegalArgumentException expected");
        } catch (MathIllegalArgumentException ex) {
            // ignored
        }
    }

[451, 'src/test/java', 'org.apache.commons.math3.linear', 'ConjugateGradientTest', 'testDimensionMismatchRightHandSide', 42, 50]

[451, 'src/test/java', 'org.apache.commons.math3.linear', 'ConjugateGradientTest', 'testDimensionMismatchSolution', 52, 60]

[451, 'src/test/java', 'org.apache.commons.math3.linear', 'SymmLQTest', 'testDimensionMismatchSolution', 234, 242]

    @Test(expected = DimensionMismatchException.class)
    public void testDimensionMismatchRightHandSide() {
        final Array2DRowRealMatrix a = new Array2DRowRealMatrix(3, 3);
        final IterativeLinearSolver solver;
        solver = new ConjugateGradient(10, 0., false);
        final ArrayRealVector b = new ArrayRealVector(2);
        final ArrayRealVector x = new ArrayRealVector(3);
        solver.solve(a, b, x);
    }

    @Test(expected = DimensionMismatchException.class)
    public void testDimensionMismatchSolution() {
        final Array2DRowRealMatrix a = new Array2DRowRealMatrix(3, 3);
        final IterativeLinearSolver solver;
        solver = new ConjugateGradient(10, 0., false);
        final ArrayRealVector b = new ArrayRealVector(3);
        final ArrayRealVector x = new ArrayRealVector(2);
        solver.solve(a, b, x);
    }

    @Test(expected = DimensionMismatchException.class)
    public void testDimensionMismatchSolution() {
        final Array2DRowRealMatrix a = new Array2DRowRealMatrix(3, 3);
        final IterativeLinearSolver solver;
        solver = new SymmLQ(10, 0., false);
        final ArrayRealVector b = new ArrayRealVector(3);
        final ArrayRealVector x = new ArrayRealVector(2);
        solver.solve(a, b, x);
    }

[478, 'src/test/java', 'org.apache.commons.math3.fraction', 'BigFractionFormatTest', 'testParseInvalid', 139, 155]

[478, 'src/test/java', 'org.apache.commons.math3.fraction', 'BigFractionFormatTest', 'testParseInvalidDenominator', 157, 173]

[478, 'src/test/java', 'org.apache.commons.math3.fraction', 'FractionFormatTest', 'testParseInvalid', 177, 193]

[478, 'src/test/java', 'org.apache.commons.math3.fraction', 'FractionFormatTest', 'testParseInvalidDenominator', 195, 211]

    @Test
    public void testParseInvalid() {
        String source = "a";
        String msg = "should not be able to parse '10 / a'.";
        try {
            properFormat.parse(source);
            Assert.fail(msg);
        } catch (MathParseException ex) {
            // success
        }
        try {
            improperFormat.parse(source);
            Assert.fail(msg);
        } catch (MathParseException ex) {
            // success
        }
    }

    @Test
    public void testParseInvalidDenominator() {
        String source = "10 / a";
        String msg = "should not be able to parse '10 / a'.";
        try {
            properFormat.parse(source);
            Assert.fail(msg);
        } catch (MathParseException ex) {
            // success
        }
        try {
            improperFormat.parse(source);
            Assert.fail(msg);
        } catch (MathParseException ex) {
            // success
        }
    }

    @Test
    public void testParseInvalid() {
        String source = "a";
        String msg = "should not be able to parse '10 / a'.";
        try {
            properFormat.parse(source);
            Assert.fail(msg);
        } catch (MathParseException ex) {
            // success
        }
        try {
            improperFormat.parse(source);
            Assert.fail(msg);
        } catch (MathParseException ex) {
            // success
        }
    }

    @Test
    public void testParseInvalidDenominator() {
        String source = "10 / a";
        String msg = "should not be able to parse '10 / a'.";
        try {
            properFormat.parse(source);
            Assert.fail(msg);
        } catch (MathParseException ex) {
            // success
        }
        try {
            improperFormat.parse(source);
            Assert.fail(msg);
        } catch (MathParseException ex) {
            // success
        }
    }

[480, 'src/test/java', 'org.apache.commons.math3.complex', 'ComplexFormatAbstractTest', 'testNan', 154, 160]

[480, 'src/test/java', 'org.apache.commons.math3.complex', 'ComplexFormatAbstractTest', 'testPositiveInfinity', 162, 168]

[480, 'src/test/java', 'org.apache.commons.math3.complex', 'ComplexFormatAbstractTest', 'testNegativeInfinity', 170, 176]

[480, 'src/test/java', 'org.apache.commons.math3.complex', 'ComplexFormatAbstractTest', 'testParseNan', 250, 256]

[480, 'src/test/java', 'org.apache.commons.math3.complex', 'ComplexFormatAbstractTest', 'testParsePositiveInfinity', 258, 264]

[480, 'src/test/java', 'org.apache.commons.math3.complex', 'ComplexFormatAbstractTest', 'testPaseNegativeInfinity', 266, 272]

    @Test
    public void testNan() {
        Complex c = new Complex(Double.NaN, Double.NaN);
        String expected = "(NaN) + (NaN)i";
        String actual = complexFormat.format(c);
        Assert.assertEquals(expected, actual);
    }

    @Test
    public void testPositiveInfinity() {
        Complex c = new Complex(Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY);
        String expected = "(Infinity) + (Infinity)i";
        String actual = complexFormat.format(c);
        Assert.assertEquals(expected, actual);
    }

    @Test
    public void testNegativeInfinity() {
        Complex c = new Complex(Double.NEGATIVE_INFINITY, Double.NEGATIVE_INFINITY);
        String expected = "(-Infinity) - (Infinity)i";
        String actual = complexFormat.format(c);
        Assert.assertEquals(expected, actual);
    }

    @Test
    public void testParseNan() {
        String source = "(NaN) + (NaN)i";
        Complex expected = new Complex(Double.NaN, Double.NaN);
        Complex actual = complexFormat.parse(source);
        Assert.assertEquals(expected, actual);
    }

    @Test
    public void testParsePositiveInfinity() {
        String source = "(Infinity) + (Infinity)i";
        Complex expected = new Complex(Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY);
        Complex actual = complexFormat.parse(source);
        Assert.assertEquals(expected, actual);
    }

    @Test
    public void testPaseNegativeInfinity() {
        String source = "(-Infinity) - (Infinity)i";
        Complex expected = new Complex(Double.NEGATIVE_INFINITY, Double.NEGATIVE_INFINITY);
        Complex actual = complexFormat.parse(source);
        Assert.assertEquals(expected, actual);
    }

[484, 'src/test/java', 'org.apache.commons.math3.geometry.euclidean.oned', 'Vector1DFormatAbstractTest', 'testParseSimpleWithDecimals', 150, 158]

[484, 'src/test/java', 'org.apache.commons.math3.geometry.euclidean.oned', 'Vector1DFormatAbstractTest', 'testParseSimpleWithDecimalsTrunc', 160, 168]

[484, 'src/test/java', 'org.apache.commons.math3.geometry.euclidean.oned', 'Vector1DFormatAbstractTest', 'testParseNegativeY', 180, 188]

[484, 'src/test/java', 'org.apache.commons.math3.geometry.euclidean.oned', 'Vector1DFormatAbstractTest', 'testParseNegativeZ', 190, 198]

[484, 'src/test/java', 'org.apache.commons.math3.geometry.euclidean.oned', 'Vector1DFormatAbstractTest', 'testParseZeroX', 210, 218]

[484, 'src/test/java', 'org.apache.commons.math3.geometry.euclidean.oned', 'Vector1DFormatAbstractTest', 'testParseNonDefaultSetting', 220, 228]

    @Test
    public void testParseSimpleWithDecimals() throws MathParseException {
        String source =
            "{1" + getDecimalCharacter() +
            "23}";
        Vector1D expected = new Vector1D(1.23);
        Vector1D actual = vector1DFormat.parse(source);
        Assert.assertEquals(expected, actual);
    }

    @Test
    public void testParseSimpleWithDecimalsTrunc() throws MathParseException {
        String source =
            "{1" + getDecimalCharacter() +
            "2323}";
        Vector1D expected = new Vector1D(1.2323);
        Vector1D actual = vector1DFormat.parse(source);
        Assert.assertEquals(expected, actual);
    }

    @Test
    public void testParseNegativeY() throws MathParseException {
        String source =
            "{1" + getDecimalCharacter() +
            "2323}";
        Vector1D expected = new Vector1D(1.2323);
        Vector1D actual = vector1DFormat.parse(source);
        Assert.assertEquals(expected, actual);
    }

    @Test
    public void testParseNegativeZ() throws MathParseException {
        String source =
            "{1" + getDecimalCharacter() +
            "2323}";
        Vector1D expected = new Vector1D(1.2323);
        Vector1D actual = vector1DFormat.parse(source);
        Assert.assertEquals(expected, actual);
    }

    @Test
    public void testParseZeroX() throws MathParseException {
        String source =
            "{0" + getDecimalCharacter() +
            "0}";
        Vector1D expected = new Vector1D(0.0);
        Vector1D actual = vector1DFormat.parse(source);
        Assert.assertEquals(expected, actual);
    }

    @Test
    public void testParseNonDefaultSetting() throws MathParseException {
        String source =
            "[1" + getDecimalCharacter() +
            "2323]";
        Vector1D expected = new Vector1D(1.2323);
        Vector1D actual = vector1DFormatSquare.parse(source);
        Assert.assertEquals(expected, actual);
    }

[496, 'src/test/java', 'org.apache.commons.math3.geometry.euclidean.threed', 'Vector3DFormatAbstractTest', 'testParseNegativeX', 208, 218]

[496, 'src/test/java', 'org.apache.commons.math3.geometry.euclidean.threed', 'Vector3DFormatAbstractTest', 'testParseNegativeY', 220, 230]

[496, 'src/test/java', 'org.apache.commons.math3.geometry.euclidean.threed', 'Vector3DFormatAbstractTest', 'testParseNegativeZ', 232, 242]

[496, 'src/test/java', 'org.apache.commons.math3.geometry.euclidean.threed', 'Vector3DFormatAbstractTest', 'testParseZeroX', 256, 266]

    @Test
    public void testParseNegativeX() throws MathParseException {
        String source =
            "{-1" + getDecimalCharacter() +
            "2323; 1" + getDecimalCharacter() +
            "4343; 1" + getDecimalCharacter() +
            "6333}";
        Vector3D expected = new Vector3D(-1.2323, 1.4343, 1.6333);
        Vector3D actual = vector3DFormat.parse(source);
        Assert.assertEquals(expected, actual);
    }

    @Test
    public void testParseNegativeY() throws MathParseException {
        String source =
            "{1" + getDecimalCharacter() +
            "2323; -1" + getDecimalCharacter() +
            "4343; 1" + getDecimalCharacter() +
            "6333}";
        Vector3D expected = new Vector3D(1.2323, -1.4343, 1.6333);
        Vector3D actual = vector3DFormat.parse(source);
        Assert.assertEquals(expected, actual);
    }

    @Test
    public void testParseNegativeZ() throws MathParseException {
        String source =
            "{1" + getDecimalCharacter() +
            "2323; 1" + getDecimalCharacter() +
            "4343; -1" + getDecimalCharacter() +
            "6333}";
        Vector3D expected = new Vector3D(1.2323, 1.4343, -1.6333);
        Vector3D actual = vector3DFormat.parse(source);
        Assert.assertEquals(expected, actual);
    }

    @Test
    public void testParseZeroX() throws MathParseException {
        String source =
            "{0" + getDecimalCharacter() +
            "0; -1" + getDecimalCharacter() +
            "4343; 1" + getDecimalCharacter() +
            "6333}";
        Vector3D expected = new Vector3D(0.0, -1.4343, 1.6333);
        Vector3D actual = vector3DFormat.parse(source);
        Assert.assertEquals(expected, actual);
    }

[504, 'src/test/java', 'org.apache.commons.math3.distribution', 'ExponentialDistributionTest', 'testInverseCumulativeProbabilityExtremes', 79, 84]

[504, 'src/test/java', 'org.apache.commons.math3.distribution', 'FDistributionTest', 'testInverseCumulativeProbabilityExtremes', 76, 81]

[504, 'src/test/java', 'org.apache.commons.math3.distribution', 'GammaDistributionTest', 'testInverseCumulativeProbabilityExtremes', 158, 163]

[504, 'src/test/java', 'org.apache.commons.math3.distribution', 'LogNormalDistributionTest', 'testInverseCumulativeProbabilityExtremes', 157, 163]

    @Test
    public void testInverseCumulativeProbabilityExtremes() {
         setInverseCumulativeTestPoints(new double[] {0, 1});
         setInverseCumulativeTestValues(new double[] {0, Double.POSITIVE_INFINITY});
         verifyInverseCumulativeProbabilities();
    }

    @Test
    public void testInverseCumulativeProbabilityExtremes() {
        setInverseCumulativeTestPoints(new double[] {0, 1});
        setInverseCumulativeTestValues(new double[] {0, Double.POSITIVE_INFINITY});
        verifyInverseCumulativeProbabilities();
    }

    @Test
    public void testInverseCumulativeProbabilityExtremes() {
        setInverseCumulativeTestPoints(new double[] {0, 1});
        setInverseCumulativeTestValues(new double[] {0, Double.POSITIVE_INFINITY});
        verifyInverseCumulativeProbabilities();
    }

    @Test
    public void testInverseCumulativeProbabilityExtremes() {
        setInverseCumulativeTestPoints(new double[] {0, 1});
        setInverseCumulativeTestValues(
                new double[] {0, Double.POSITIVE_INFINITY});
        verifyInverseCumulativeProbabilities();
    }

[505, 'src/test/java', 'org.apache.commons.math3.analysis.function', 'GaussianTest', 'testParametricUsage2', 101, 105]

[505, 'src/test/java', 'org.apache.commons.math3.analysis.function', 'GaussianTest', 'testParametricUsage5', 119, 123]

[505, 'src/test/java', 'org.apache.commons.math3.analysis.function', 'HarmonicOscillatorTest', 'testParametricUsage2', 91, 95]

[505, 'src/test/java', 'org.apache.commons.math3.analysis.function', 'HarmonicOscillatorTest', 'testParametricUsage4', 103, 107]

[505, 'src/test/java', 'org.apache.commons.math3.analysis.function', 'LogisticTest', 'testParametricUsage2', 112, 116]

[505, 'src/test/java', 'org.apache.commons.math3.analysis.function', 'LogisticTest', 'testParametricUsage4', 124, 128]

[505, 'src/test/java', 'org.apache.commons.math3.analysis.function', 'LogitTest', 'testParametricUsage2', 114, 118]

[505, 'src/test/java', 'org.apache.commons.math3.analysis.function', 'LogitTest', 'testParametricUsage4', 126, 130]

[505, 'src/test/java', 'org.apache.commons.math3.analysis.function', 'SigmoidTest', 'testParametricUsage2', 82, 86]

[505, 'src/test/java', 'org.apache.commons.math3.analysis.function', 'SigmoidTest', 'testParametricUsage4', 94, 98]

    @Test(expected=DimensionMismatchException.class)
    public void testParametricUsage2() {
        final Gaussian.Parametric g = new Gaussian.Parametric();
        g.value(0, new double[] {0});
    }

    @Test(expected=DimensionMismatchException.class)
    public void testParametricUsage5() {
        final Gaussian.Parametric g = new Gaussian.Parametric();
        g.gradient(0, new double[] {0});
    }

    @Test(expected=DimensionMismatchException.class)
    public void testParametricUsage2() {
        final HarmonicOscillator.Parametric g = new HarmonicOscillator.Parametric();
        g.value(0, new double[] {0});
    }

    @Test(expected=DimensionMismatchException.class)
    public void testParametricUsage4() {
        final HarmonicOscillator.Parametric g = new HarmonicOscillator.Parametric();
        g.gradient(0, new double[] {0});
    }

    @Test(expected=DimensionMismatchException.class)
    public void testParametricUsage2() {
        final Logistic.Parametric g = new Logistic.Parametric();
        g.value(0, new double[] {0});
    }

    @Test(expected=DimensionMismatchException.class)
    public void testParametricUsage4() {
        final Logistic.Parametric g = new Logistic.Parametric();
        g.gradient(0, new double[] {0});
    }

    @Test(expected=DimensionMismatchException.class)
    public void testParametricUsage2() {
        final Logit.Parametric g = new Logit.Parametric();
        g.value(0, new double[] {0});
    }

    @Test(expected=DimensionMismatchException.class)
    public void testParametricUsage4() {
        final Logit.Parametric g = new Logit.Parametric();
        g.gradient(0, new double[] {0});
    }

    @Test(expected=DimensionMismatchException.class)
    public void testParametricUsage2() {
        final Sigmoid.Parametric g = new Sigmoid.Parametric();
        g.value(0, new double[] {0});
    }

    @Test(expected=DimensionMismatchException.class)
    public void testParametricUsage4() {
        final Sigmoid.Parametric g = new Sigmoid.Parametric();
        g.gradient(0, new double[] {0});
    }

[507, 'src/test/java', 'org.apache.commons.math3.linear', 'RealVectorAbstractTest', 'testProjectionNullVector', 978, 981]

[507, 'src/test/java', 'org.apache.commons.math3.linear', 'RealVectorAbstractTest', 'testGetDistanceDimensionMismatch', 629, 632]

[507, 'src/test/java', 'org.apache.commons.math3.linear', 'RealVectorAbstractTest', 'testGetL1DistanceDimensionMismatch', 678, 681]

[507, 'src/test/java', 'org.apache.commons.math3.linear', 'RealVectorAbstractTest', 'testGetLInfDistanceDimensionMismatch', 727, 730]

    @Test(expected = MathArithmeticException.class)
    public void testProjectionNullVector() {
        create(new double[4]).projection(create(new double[4]));
    }

    @Test(expected = DimensionMismatchException.class)
    public void testGetDistanceDimensionMismatch() {
        create(new double[4]).getDistance(createAlien(new double[5]));
    }

    @Test(expected = DimensionMismatchException.class)
    public void testGetL1DistanceDimensionMismatch() {
        create(new double[4]).getL1Distance(createAlien(new double[5]));
    }

    @Test(expected = DimensionMismatchException.class)
    public void testGetLInfDistanceDimensionMismatch() {
        create(new double[4]).getLInfDistance(createAlien(new double[5]));
    }

[520, 'src/test/java', 'org.apache.commons.math3.linear', 'BlockFieldMatrixTest', 'testGetSubMatrix', 568, 586]

[520, 'src/test/java', 'org.apache.commons.math3.linear', 'BlockFieldMatrixTest', 'testCopySubMatrix', 666, 685]

[520, 'src/test/java', 'org.apache.commons.math3.linear', 'FieldMatrixImplTest', 'testGetSubMatrix', 428, 446]

[520, 'src/test/java', 'org.apache.commons.math3.linear', 'FieldMatrixImplTest', 'testCopySubMatrix', 506, 525]

    @Test
    public void testGetSubMatrix() {
        FieldMatrix<Fraction> m = new BlockFieldMatrix<Fraction>(subTestData);
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
        FieldMatrix<Fraction> m = new BlockFieldMatrix<Fraction>(subTestData);
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
        checkCopy(m, null, new int[] {}, new int[] { 0 });
        checkCopy(m, null, new int[] { 0 }, new int[] { 4 });
    }

    @Test
    public void testGetSubMatrix() {
        FieldMatrix<Fraction> m = new Array2DRowFieldMatrix<Fraction>(subTestData);
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
        FieldMatrix<Fraction> m = new Array2DRowFieldMatrix<Fraction>(subTestData);
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

[612, 'src/test/java', 'org.apache.commons.math3.stat.interval', 'IntervalUtilsTest', 'testAgrestiCoull', 35, 39]

[612, 'src/test/java', 'org.apache.commons.math3.stat.interval', 'IntervalUtilsTest', 'testClopperPearson', 41, 45]

[612, 'src/test/java', 'org.apache.commons.math3.stat.interval', 'IntervalUtilsTest', 'testNormalApproximation', 47, 51]

[612, 'src/test/java', 'org.apache.commons.math3.stat.interval', 'IntervalUtilsTest', 'testWilsonScore', 53, 57]

    @Test
    public void testAgrestiCoull() {
        checkConfidenceIntervals(new AgrestiCoullInterval().createInterval(trials, successes, confidenceLevel),
                                 IntervalUtils.getAgrestiCoullInterval(trials, successes, confidenceLevel));
    }

    @Test
    public void testClopperPearson() {
        checkConfidenceIntervals(new ClopperPearsonInterval().createInterval(trials, successes, confidenceLevel),
                                 IntervalUtils.getClopperPearsonInterval(trials, successes, confidenceLevel));
    }

    @Test
    public void testNormalApproximation() {
        checkConfidenceIntervals(new NormalApproximationInterval().createInterval(trials, successes, confidenceLevel),
                                 IntervalUtils.getNormalApproximationInterval(trials, successes, confidenceLevel));
    }

    @Test
    public void testWilsonScore() {
        checkConfidenceIntervals(new WilsonScoreInterval().createInterval(trials, successes, confidenceLevel),
                                 IntervalUtils.getWilsonScoreInterval(trials, successes, confidenceLevel));
    }

[662, 'src/test/java', 'org.apache.commons.math3.analysis.differentiation', 'SparseGradientTest', 'testVarAddition', 55, 69]

[662, 'src/test/java', 'org.apache.commons.math3.analysis.differentiation', 'SparseGradientTest', 'testSubtraction', 71, 85]

[662, 'src/test/java', 'org.apache.commons.math3.analysis.differentiation', 'SparseGradientTest', 'testDivision', 87, 100]

    @Test
    public void testVarAddition() {
        final double v1 = 1.0;
        final double v2 = 2.0;
        final int id1 = -1;
        final int id2 = 3;
        final SparseGradient var1 = SparseGradient.createVariable(id1, v1);
        final SparseGradient var2 = SparseGradient.createVariable(id2, v2);
        final SparseGradient sum = var1.add(var2);

        Assert.assertEquals(v1 + v2, sum.getValue(), 1.0e-15); // returns the value
        Assert.assertEquals(2, sum.numVars());
        Assert.assertEquals(1.0, sum.getDerivative(id1), 1.0e-15);
        Assert.assertEquals(1.0, sum.getDerivative(id2), 1.0e-15);
    }

    @Test
    public void testSubtraction() {
        final double v1 = 1.0;
        final double v2 = 2.0;
        final int id1 = -1;
        final int id2 = 3;
        final SparseGradient var1 = SparseGradient.createVariable(id1, v1);
        final SparseGradient var2 = SparseGradient.createVariable(id2, v2);
        final SparseGradient sum = var1.subtract(var2);

        Assert.assertEquals(v1 - v2, sum.getValue(), 1.0e-15); // returns the value
        Assert.assertEquals(2, sum.numVars());
        Assert.assertEquals(1.0, sum.getDerivative(id1), 1.0e-15);
        Assert.assertEquals(-1.0, sum.getDerivative(id2), 1.0e-15);
    }

    @Test
    public void testDivision() {
        final double v1 = 1.0;
        final double v2 = 2.0;
        final int id1 = -1;
        final int id2 = 3;
        final SparseGradient var1 = SparseGradient.createVariable(id1, v1);
        final SparseGradient var2 = SparseGradient.createVariable(id2, v2);
        final SparseGradient out = var1.divide(var2);
        Assert.assertEquals(v1 / v2, out.getValue(), 1.0e-15); // returns the value
        Assert.assertEquals(2, out.numVars());
        Assert.assertEquals(1 / v2, out.getDerivative(id1), 1.0e-15);
        Assert.assertEquals(-1 / (v2 * v2), out.getDerivative(id2), 1.0e-15);
    }

[715, 'src/test/java', 'org.apache.commons.math3.ode.nonstiff', 'ClassicalRungeKuttaIntegratorTest', 'testSmallStep', 184, 202]

[715, 'src/test/java', 'org.apache.commons.math3.ode.nonstiff', 'EulerIntegratorTest', 'testSmallStep', 103, 123]

[715, 'src/test/java', 'org.apache.commons.math3.ode.nonstiff', 'GillIntegratorTest', 'testSmallStep', 103, 122]

[715, 'src/test/java', 'org.apache.commons.math3.ode.nonstiff', 'LutherIntegratorTest', 'testSmallStep', 184, 202]

[715, 'src/test/java', 'org.apache.commons.math3.ode.nonstiff', 'MidpointIntegratorTest', 'testSmallStep', 103, 123]

[715, 'src/test/java', 'org.apache.commons.math3.ode.nonstiff', 'ThreeEighthesIntegratorTest', 'testSmallStep', 103, 122]

  @Test
  public void testSmallStep()
      throws DimensionMismatchException, NumberIsTooSmallException,
             MaxCountExceededException, NoBracketingException {

    TestProblem1 pb = new TestProblem1();
    double step = (pb.getFinalTime() - pb.getInitialTime()) * 0.001;

    FirstOrderIntegrator integ = new ClassicalRungeKuttaIntegrator(step);
    TestProblemHandler handler = new TestProblemHandler(pb, integ);
    integ.addStepHandler(handler);
    integ.integrate(pb, pb.getInitialTime(), pb.getInitialState(),
                    pb.getFinalTime(), new double[pb.getDimension()]);

    Assert.assertTrue(handler.getLastError() < 2.0e-13);
    Assert.assertTrue(handler.getMaximalValueError() < 4.0e-12);
    Assert.assertEquals(0, handler.getMaximalTimeError(), 1.0e-12);
    Assert.assertEquals("classical Runge-Kutta", integ.getName());
  }

  @Test
  public void testSmallStep()
      throws DimensionMismatchException, NumberIsTooSmallException,
             MaxCountExceededException, NoBracketingException {

    TestProblem1 pb  = new TestProblem1();
    double step = (pb.getFinalTime() - pb.getInitialTime()) * 0.001;

    FirstOrderIntegrator integ = new EulerIntegrator(step);
    TestProblemHandler handler = new TestProblemHandler(pb, integ);
    integ.addStepHandler(handler);
    integ.integrate(pb,
                    pb.getInitialTime(), pb.getInitialState(),
                    pb.getFinalTime(), new double[pb.getDimension()]);

   Assert.assertTrue(handler.getLastError() < 2.0e-4);
   Assert.assertTrue(handler.getMaximalValueError() < 1.0e-3);
   Assert.assertEquals(0, handler.getMaximalTimeError(), 1.0e-12);
   Assert.assertEquals("Euler", integ.getName());

  }

  @Test
  public void testSmallStep()
      throws DimensionMismatchException, NumberIsTooSmallException,
             MaxCountExceededException, NoBracketingException {

    TestProblem1 pb = new TestProblem1();
    double step = (pb.getFinalTime() - pb.getInitialTime()) * 0.001;

    FirstOrderIntegrator integ = new GillIntegrator(step);
    TestProblemHandler handler = new TestProblemHandler(pb, integ);
    integ.addStepHandler(handler);
    integ.integrate(pb, pb.getInitialTime(), pb.getInitialState(),
                    pb.getFinalTime(), new double[pb.getDimension()]);

    Assert.assertTrue(handler.getLastError() < 2.0e-13);
    Assert.assertTrue(handler.getMaximalValueError() < 4.0e-12);
    Assert.assertEquals(0, handler.getMaximalTimeError(), 1.0e-12);
    Assert.assertEquals("Gill", integ.getName());

  }

    @Test
    public void testSmallStep()
            throws DimensionMismatchException, NumberIsTooSmallException,
            MaxCountExceededException, NoBracketingException {

        TestProblem1 pb = new TestProblem1();
        double step = (pb.getFinalTime() - pb.getInitialTime()) * 0.001;

        FirstOrderIntegrator integ = new LutherIntegrator(step);
        TestProblemHandler handler = new TestProblemHandler(pb, integ);
        integ.addStepHandler(handler);
        integ.integrate(pb, pb.getInitialTime(), pb.getInitialState(),
                        pb.getFinalTime(), new double[pb.getDimension()]);

        Assert.assertTrue(handler.getLastError() < 9.0e-17);
        Assert.assertTrue(handler.getMaximalValueError() < 4.0e-15);
        Assert.assertEquals(0, handler.getMaximalTimeError(), 1.0e-12);
        Assert.assertEquals("Luther", integ.getName());
    }

  @Test
  public void testSmallStep()
      throws DimensionMismatchException, NumberIsTooSmallException,
             MaxCountExceededException, NoBracketingException {

    TestProblem1 pb  = new TestProblem1();
    double step = (pb.getFinalTime() - pb.getInitialTime()) * 0.001;

    FirstOrderIntegrator integ = new MidpointIntegrator(step);
    TestProblemHandler handler = new TestProblemHandler(pb, integ);
    integ.addStepHandler(handler);
    integ.integrate(pb,
                    pb.getInitialTime(), pb.getInitialState(),
                    pb.getFinalTime(), new double[pb.getDimension()]);

    Assert.assertTrue(handler.getLastError() < 2.0e-7);
    Assert.assertTrue(handler.getMaximalValueError() < 1.0e-6);
    Assert.assertEquals(0, handler.getMaximalTimeError(), 1.0e-12);
    Assert.assertEquals("midpoint", integ.getName());

  }

 @Test
 public void testSmallStep()
     throws DimensionMismatchException, NumberIsTooSmallException,
            MaxCountExceededException, NoBracketingException {

    TestProblem1 pb = new TestProblem1();
    double step = (pb.getFinalTime() - pb.getInitialTime()) * 0.001;

    FirstOrderIntegrator integ = new ThreeEighthesIntegrator(step);
    TestProblemHandler handler = new TestProblemHandler(pb, integ);
    integ.addStepHandler(handler);
    integ.integrate(pb, pb.getInitialTime(), pb.getInitialState(),
                    pb.getFinalTime(), new double[pb.getDimension()]);

    Assert.assertTrue(handler.getLastError() < 2.0e-13);
    Assert.assertTrue(handler.getMaximalValueError() < 4.0e-12);
    Assert.assertEquals(0, handler.getMaximalTimeError(), 1.0e-12);
    Assert.assertEquals("3/8", integ.getName());

  }

[720, 'src/test/java', 'org.apache.commons.math3.linear', 'ConjugateGradientTest', 'testNonSquarePreconditioner', 219, 243]

[720, 'src/test/java', 'org.apache.commons.math3.linear', 'ConjugateGradientTest', 'testMismatchedOperatorDimensions', 245, 269]

[720, 'src/test/java', 'org.apache.commons.math3.linear', 'SymmLQTest', 'testNonSquarePreconditioner', 320, 344]

[720, 'src/test/java', 'org.apache.commons.math3.linear', 'SymmLQTest', 'testMismatchedOperatorDimensions', 346, 370]

    @Test(expected = NonSquareOperatorException.class)
    public void testNonSquarePreconditioner() {
        final Array2DRowRealMatrix a = new Array2DRowRealMatrix(2, 2);
        final RealLinearOperator m = new RealLinearOperator() {

            @Override
            public RealVector operate(final RealVector x) {
                throw new UnsupportedOperationException();
            }

            @Override
            public int getRowDimension() {
                return 2;
            }

            @Override
            public int getColumnDimension() {
                return 3;
            }
        };
        final PreconditionedIterativeLinearSolver solver;
        solver = new ConjugateGradient(10, 0d, false);
        final ArrayRealVector b = new ArrayRealVector(a.getRowDimension());
        solver.solve(a, m, b);
    }

    @Test(expected = DimensionMismatchException.class)
    public void testMismatchedOperatorDimensions() {
        final Array2DRowRealMatrix a = new Array2DRowRealMatrix(2, 2);
        final RealLinearOperator m = new RealLinearOperator() {

            @Override
            public RealVector operate(final RealVector x) {
                throw new UnsupportedOperationException();
            }

            @Override
            public int getRowDimension() {
                return 3;
            }

            @Override
            public int getColumnDimension() {
                return 3;
            }
        };
        final PreconditionedIterativeLinearSolver solver;
        solver = new ConjugateGradient(10, 0d, false);
        final ArrayRealVector b = new ArrayRealVector(a.getRowDimension());
        solver.solve(a, m, b);
    }

    @Test(expected = NonSquareOperatorException.class)
    public void testNonSquarePreconditioner() {
        final Array2DRowRealMatrix a = new Array2DRowRealMatrix(2, 2);
        final RealLinearOperator m = new RealLinearOperator() {

            @Override
            public RealVector operate(final RealVector x) {
                throw new UnsupportedOperationException();
            }

            @Override
            public int getRowDimension() {
                return 2;
            }

            @Override
            public int getColumnDimension() {
                return 3;
            }
        };
        final PreconditionedIterativeLinearSolver solver;
        solver = new SymmLQ(10, 0., false);
        final ArrayRealVector b = new ArrayRealVector(a.getRowDimension());
        solver.solve(a, m, b);
    }

    @Test(expected = DimensionMismatchException.class)
    public void testMismatchedOperatorDimensions() {
        final Array2DRowRealMatrix a = new Array2DRowRealMatrix(2, 2);
        final RealLinearOperator m = new RealLinearOperator() {

            @Override
            public RealVector operate(final RealVector x) {
                throw new UnsupportedOperationException();
            }

            @Override
            public int getRowDimension() {
                return 3;
            }

            @Override
            public int getColumnDimension() {
                return 3;
            }
        };
        final PreconditionedIterativeLinearSolver solver;
        solver = new SymmLQ(10, 0d, false);
        final ArrayRealVector b = new ArrayRealVector(a.getRowDimension());
        solver.solve(a, m, b);
    }

[722, 'src/test/java', 'org.apache.commons.math3.complex', 'ComplexTest', 'testAcosInf', 665, 675]

[722, 'src/test/java', 'org.apache.commons.math3.complex', 'ComplexTest', 'testAsinInf', 694, 704]

[722, 'src/test/java', 'org.apache.commons.math3.complex', 'ComplexTest', 'testAtanInf', 714, 724]

    @Test
    public void testAcosInf() {
        TestUtils.assertSame(Complex.NaN, oneInf.acos());
        TestUtils.assertSame(Complex.NaN, oneNegInf.acos());
        TestUtils.assertSame(Complex.NaN, infOne.acos());
        TestUtils.assertSame(Complex.NaN, negInfOne.acos());
        TestUtils.assertSame(Complex.NaN, infInf.acos());
        TestUtils.assertSame(Complex.NaN, infNegInf.acos());
        TestUtils.assertSame(Complex.NaN, negInfInf.acos());
        TestUtils.assertSame(Complex.NaN, negInfNegInf.acos());
    }

    @Test
    public void testAsinInf() {
        TestUtils.assertSame(Complex.NaN, oneInf.asin());
        TestUtils.assertSame(Complex.NaN, oneNegInf.asin());
        TestUtils.assertSame(Complex.NaN, infOne.asin());
        TestUtils.assertSame(Complex.NaN, negInfOne.asin());
        TestUtils.assertSame(Complex.NaN, infInf.asin());
        TestUtils.assertSame(Complex.NaN, infNegInf.asin());
        TestUtils.assertSame(Complex.NaN, negInfInf.asin());
        TestUtils.assertSame(Complex.NaN, negInfNegInf.asin());
    }

    @Test
    public void testAtanInf() {
        TestUtils.assertSame(Complex.NaN, oneInf.atan());
        TestUtils.assertSame(Complex.NaN, oneNegInf.atan());
        TestUtils.assertSame(Complex.NaN, infOne.atan());
        TestUtils.assertSame(Complex.NaN, negInfOne.atan());
        TestUtils.assertSame(Complex.NaN, infInf.atan());
        TestUtils.assertSame(Complex.NaN, infNegInf.atan());
        TestUtils.assertSame(Complex.NaN, negInfInf.atan());
        TestUtils.assertSame(Complex.NaN, negInfNegInf.atan());
    }

[728, 'src/test/java', 'org.apache.commons.math3.linear', 'RealMatrixFormatAbstractTest', 'testNan', 158, 164]

[728, 'src/test/java', 'org.apache.commons.math3.linear', 'RealMatrixFormatAbstractTest', 'testPositiveInfinity', 166, 173]

[728, 'src/test/java', 'org.apache.commons.math3.linear', 'RealMatrixFormatAbstractTest', 'tesNegativeInfinity', 175, 182]

[728, 'src/test/java', 'org.apache.commons.math3.linear', 'RealMatrixFormatAbstractTest', 'testParsePositiveInfinity', 298, 305]

[728, 'src/test/java', 'org.apache.commons.math3.linear', 'RealMatrixFormatAbstractTest', 'testParseNegativeInfinity', 307, 314]

    @Test
    public void testNan() {
        RealMatrix m = MatrixUtils.createRealMatrix(new double[][] {{Double.NaN, Double.NaN, Double.NaN}});
        String expected = "{{(NaN),(NaN),(NaN)}}";
        String actual = realMatrixFormat.format(m);
        Assert.assertEquals(expected, actual);
    }

    @Test
    public void testPositiveInfinity() {
        RealMatrix m = MatrixUtils.createRealMatrix(
                new double[][] {{Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY}});
        String expected = "{{(Infinity),(Infinity),(Infinity)}}";
        String actual = realMatrixFormat.format(m);
        Assert.assertEquals(expected, actual);
    }

    @Test
    public void tesNegativeInfinity() {
        RealMatrix m = MatrixUtils.createRealMatrix(
                new double[][] {{Double.NEGATIVE_INFINITY, Double.NEGATIVE_INFINITY, Double.NEGATIVE_INFINITY}});
        String expected = "{{(-Infinity),(-Infinity),(-Infinity)}}";
        String actual = realMatrixFormat.format(m);
        Assert.assertEquals(expected, actual);
    }

    @Test
    public void testParsePositiveInfinity() {
        String source = "{{(Infinity), (Infinity), (Infinity)}}";
        RealMatrix actual = realMatrixFormat.parse(source);
        RealMatrix expected = MatrixUtils.createRealMatrix(
                new double[][] {{Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY}});
        Assert.assertEquals(expected, actual);
    }

    @Test
    public void testParseNegativeInfinity() {
        String source = "{{(-Infinity), (-Infinity), (-Infinity)}}";
        RealMatrix actual = realMatrixFormat.parse(source);
        RealMatrix expected = MatrixUtils.createRealMatrix(
                new double[][] {{Double.NEGATIVE_INFINITY, Double.NEGATIVE_INFINITY, Double.NEGATIVE_INFINITY}});
        Assert.assertEquals(expected, actual);
    }

[759, 'src/test/java', 'org.apache.commons.math3.linear', 'RealVectorFormatAbstractTest', 'testForgottenPrefix', 332, 338]

[759, 'src/test/java', 'org.apache.commons.math3.linear', 'RealVectorFormatAbstractTest', 'testForgottenSeparator', 340, 346]

[759, 'src/test/java', 'org.apache.commons.math3.linear', 'RealVectorFormatAbstractTest', 'testForgottenSuffix', 348, 354]

    @Test
    public void testForgottenPrefix() {
        ParsePosition pos = new ParsePosition(0);
        final String source = "1; 1; 1}";
        Assert.assertNull("Should not parse <"+source+">",new RealVectorFormat().parse(source, pos));
        Assert.assertEquals(0, pos.getErrorIndex());
    }

    @Test
    public void testForgottenSeparator() {
        ParsePosition pos = new ParsePosition(0);
        final String source = "{1; 1 1}";
        Assert.assertNull("Should not parse <"+source+">",new RealVectorFormat().parse(source, pos));
        Assert.assertEquals(6, pos.getErrorIndex());
    }

    @Test
    public void testForgottenSuffix() {
        ParsePosition pos = new ParsePosition(0);
        final String source = "{1; 1; 1 ";
        Assert.assertNull("Should not parse <"+source+">",new RealVectorFormat().parse(source, pos));
        Assert.assertEquals(8, pos.getErrorIndex());
    }

[783, 'src/test/java', 'org.apache.commons.math3.linear', 'Array2DRowRealMatrixTest', 'testSetRowMatrix', 637, 656]

[783, 'src/test/java', 'org.apache.commons.math3.linear', 'Array2DRowRealMatrixTest', 'testSetColumnMatrix', 681, 700]

[783, 'src/test/java', 'org.apache.commons.math3.linear', 'BlockRealMatrixTest', 'testSetRowMatrix', 681, 700]

[783, 'src/test/java', 'org.apache.commons.math3.linear', 'BlockRealMatrixTest', 'testSetColumnMatrix', 742, 761]

    @Test
    public void testSetRowMatrix() {
        RealMatrix m = new Array2DRowRealMatrix(subTestData);
        RealMatrix mRow3 = new Array2DRowRealMatrix(subRow3);
        Assert.assertNotSame(mRow3, m.getRowMatrix(0));
        m.setRowMatrix(0, mRow3);
        Assert.assertEquals(mRow3, m.getRowMatrix(0));
        try {
            m.setRowMatrix(-1, mRow3);
            Assert.fail("Expecting OutOfRangeException");
        } catch (OutOfRangeException ex) {
            // expected
        }
        try {
            m.setRowMatrix(0, m);
            Assert.fail("Expecting MatrixDimensionMismatchException");
        } catch (MatrixDimensionMismatchException ex) {
            // expected
        }
    }

    @Test
    public void testSetColumnMatrix() {
        RealMatrix m = new Array2DRowRealMatrix(subTestData);
        RealMatrix mColumn3 = new Array2DRowRealMatrix(subColumn3);
        Assert.assertNotSame(mColumn3, m.getColumnMatrix(1));
        m.setColumnMatrix(1, mColumn3);
        Assert.assertEquals(mColumn3, m.getColumnMatrix(1));
        try {
            m.setColumnMatrix(-1, mColumn3);
            Assert.fail("Expecting OutOfRangeException");
        } catch (OutOfRangeException ex) {
            // expected
        }
        try {
            m.setColumnMatrix(0, m);
            Assert.fail("Expecting MatrixDimensionMismatchException");
        } catch (MatrixDimensionMismatchException ex) {
            // expected
        }
    }

    @Test
    public void testSetRowMatrix() {
        RealMatrix m = new BlockRealMatrix(subTestData);
        RealMatrix mRow3 = new BlockRealMatrix(subRow3);
        Assert.assertNotSame(mRow3, m.getRowMatrix(0));
        m.setRowMatrix(0, mRow3);
        Assert.assertEquals(mRow3, m.getRowMatrix(0));
        try {
            m.setRowMatrix(-1, mRow3);
            Assert.fail("Expecting OutOfRangeException");
        } catch (OutOfRangeException ex) {
            // expected
        }
        try {
            m.setRowMatrix(0, m);
            Assert.fail("Expecting MatrixDimensionMismatchException");
        } catch (MatrixDimensionMismatchException ex) {
            // expected
        }
    }

    @Test
    public void testSetColumnMatrix() {
        RealMatrix m = new BlockRealMatrix(subTestData);
        RealMatrix mColumn3 = new BlockRealMatrix(subColumn3);
        Assert.assertNotSame(mColumn3, m.getColumnMatrix(1));
        m.setColumnMatrix(1, mColumn3);
        Assert.assertEquals(mColumn3, m.getColumnMatrix(1));
        try {
            m.setColumnMatrix(-1, mColumn3);
            Assert.fail("Expecting OutOfRangeException");
        } catch (OutOfRangeException ex) {
            // expected
        }
        try {
            m.setColumnMatrix(0, m);
            Assert.fail("Expecting MatrixDimensionMismatchException");
        } catch (MatrixDimensionMismatchException ex) {
            // expected
        }
    }

[794, 'src/test/java', 'org.apache.commons.math3.geometry.euclidean.twod', 'Vector2DFormatAbstractTest', 'testParseNegativeX', 198, 207]

[794, 'src/test/java', 'org.apache.commons.math3.geometry.euclidean.twod', 'Vector2DFormatAbstractTest', 'testParseNegativeY', 209, 218]

[794, 'src/test/java', 'org.apache.commons.math3.geometry.euclidean.twod', 'Vector2DFormatAbstractTest', 'testParseZeroX', 242, 251]

    @Test
    public void testParseNegativeX() throws MathParseException {
        String source =
            "{-1" + getDecimalCharacter() +
            "2323; 1" + getDecimalCharacter() +
            "4343}";
        Vector2D expected = new Vector2D(-1.2323, 1.4343);
        Vector2D actual = vector2DFormat.parse(source);
        Assert.assertEquals(expected, actual);
    }

    @Test
    public void testParseNegativeY() throws MathParseException {
        String source =
            "{1" + getDecimalCharacter() +
            "2323; -1" + getDecimalCharacter() +
            "4343}";
        Vector2D expected = new Vector2D(1.2323, -1.4343);
        Vector2D actual = vector2DFormat.parse(source);
        Assert.assertEquals(expected, actual);
    }

    @Test
    public void testParseZeroX() throws MathParseException {
        String source =
            "{0" + getDecimalCharacter() +
            "0; -1" + getDecimalCharacter() +
            "4343}";
        Vector2D expected = new Vector2D(0.0, -1.4343);
        Vector2D actual = vector2DFormat.parse(source);
        Assert.assertEquals(expected, actual);
    }

[798, 'src/test/java', 'org.apache.commons.math3.analysis.integration', 'MidPointIntegratorTest', 'testParameters', 124, 149]

[798, 'src/test/java', 'org.apache.commons.math3.analysis.integration', 'RombergIntegratorTest', 'testParameters', 97, 122]

[798, 'src/test/java', 'org.apache.commons.math3.analysis.integration', 'SimpsonIntegratorTest', 'testParameters', 96, 120]

[798, 'src/test/java', 'org.apache.commons.math3.analysis.integration', 'TrapezoidIntegratorTest', 'testParameters', 97, 122]

    @Test
    public void testParameters() {
        UnivariateFunction f = new Sin();

        try {
            // bad interval
            new MidPointIntegrator().integrate(1000, f, 1, -1);
            Assert.fail("Expecting NumberIsTooLargeException - bad interval");
        } catch (NumberIsTooLargeException ex) {
            // expected
        }
        try {
            // bad iteration limits
            new MidPointIntegrator(5, 4);
            Assert.fail("Expecting NumberIsTooSmallException - bad iteration limits");
        } catch (NumberIsTooSmallException ex) {
            // expected
        }
        try {
            // bad iteration limits
            new MidPointIntegrator(10, 99);
            Assert.fail("Expecting NumberIsTooLargeException - bad iteration limits");
        } catch (NumberIsTooLargeException ex) {
            // expected
        }
    }

    @Test
    public void testParameters() {
        UnivariateFunction f = new Sin();

        try {
            // bad interval
            new RombergIntegrator().integrate(1000, f, 1, -1);
            Assert.fail("Expecting NumberIsTooLargeException - bad interval");
        } catch (NumberIsTooLargeException ex) {
            // expected
        }
        try {
            // bad iteration limits
            new RombergIntegrator(5, 4);
            Assert.fail("Expecting NumberIsTooSmallException - bad iteration limits");
        } catch (NumberIsTooSmallException ex) {
            // expected
        }
        try {
            // bad iteration limits
            new RombergIntegrator(10, 50);
            Assert.fail("Expecting NumberIsTooLargeException - bad iteration limits");
        } catch (NumberIsTooLargeException ex) {
            // expected
        }
    }

    @Test
    public void testParameters() {
        UnivariateFunction f = new Sin();
        try {
            // bad interval
            new SimpsonIntegrator().integrate(1000, f, 1, -1);
            Assert.fail("Expecting NumberIsTooLargeException - bad interval");
        } catch (NumberIsTooLargeException ex) {
            // expected
        }
        try {
            // bad iteration limits
            new SimpsonIntegrator(5, 4);
            Assert.fail("Expecting NumberIsTooSmallException - bad iteration limits");
        } catch (NumberIsTooSmallException ex) {
            // expected
        }
        try {
            // bad iteration limits
            new SimpsonIntegrator(10, 99);
            Assert.fail("Expecting NumberIsTooLargeException - bad iteration limits");
        } catch (NumberIsTooLargeException ex) {
            // expected
        }
    }

    @Test
    public void testParameters() {
        UnivariateFunction f = new Sin();

        try {
            // bad interval
            new TrapezoidIntegrator().integrate(1000, f, 1, -1);
            Assert.fail("Expecting NumberIsTooLargeException - bad interval");
        } catch (NumberIsTooLargeException ex) {
            // expected
        }
        try {
            // bad iteration limits
            new TrapezoidIntegrator(5, 4);
            Assert.fail("Expecting NumberIsTooSmallException - bad iteration limits");
        } catch (NumberIsTooSmallException ex) {
            // expected
        }
        try {
            // bad iteration limits
            new TrapezoidIntegrator(10,99);
            Assert.fail("Expecting NumberIsTooLargeException - bad iteration limits");
        } catch (NumberIsTooLargeException ex) {
            // expected
        }
    }

[862, 'src/test/java', 'org.apache.commons.math3.ode.nonstiff', 'DormandPrince54StepInterpolatorTest', 'derivativesConsistency', 43, 56]

[862, 'src/test/java', 'org.apache.commons.math3.ode.nonstiff', 'DormandPrince853StepInterpolatorTest', 'derivativesConsistency', 43, 56]

[862, 'src/test/java', 'org.apache.commons.math3.ode.nonstiff', 'HighamHall54StepInterpolatorTest', 'derivativesConsistency', 43, 56]

  @Test
  public void derivativesConsistency()
      throws DimensionMismatchException, NumberIsTooSmallException,
             MaxCountExceededException, NoBracketingException {
    TestProblem3 pb = new TestProblem3(0.1);
    double minStep = 0;
    double maxStep = pb.getFinalTime() - pb.getInitialTime();
    double scalAbsoluteTolerance = 1.0e-8;
    double scalRelativeTolerance = scalAbsoluteTolerance;
    DormandPrince54Integrator integ = new DormandPrince54Integrator(minStep, maxStep,
                                                                    scalAbsoluteTolerance,
                                                                    scalRelativeTolerance);
    StepInterpolatorTestUtils.checkDerivativesConsistency(integ, pb, 0.01, 3.6e-12);
  }

  @Test
  public void derivativesConsistency()
      throws DimensionMismatchException, NumberIsTooSmallException,
             MaxCountExceededException, NoBracketingException {
    TestProblem3 pb = new TestProblem3(0.1);
    double minStep = 0;
    double maxStep = pb.getFinalTime() - pb.getInitialTime();
    double scalAbsoluteTolerance = 1.0e-8;
    double scalRelativeTolerance = scalAbsoluteTolerance;
    DormandPrince853Integrator integ = new DormandPrince853Integrator(minStep, maxStep,
                                                                      scalAbsoluteTolerance,
                                                                      scalRelativeTolerance);
    StepInterpolatorTestUtils.checkDerivativesConsistency(integ, pb, 0.01, 1.8e-12);
  }

  @Test
  public void derivativesConsistency()
      throws DimensionMismatchException, NumberIsTooSmallException,
             MaxCountExceededException, NoBracketingException {
    TestProblem3 pb = new TestProblem3(0.1);
    double minStep = 0;
    double maxStep = pb.getFinalTime() - pb.getInitialTime();
    double scalAbsoluteTolerance = 1.0e-8;
    double scalRelativeTolerance = scalAbsoluteTolerance;
    HighamHall54Integrator integ = new HighamHall54Integrator(minStep, maxStep,
                                                              scalAbsoluteTolerance,
                                                              scalRelativeTolerance);
    StepInterpolatorTestUtils.checkDerivativesConsistency(integ, pb, 0.01, 4.8e-12);
  }

[875, 'src/test/java', 'org.apache.commons.math3.linear', 'RealVectorFormatAbstractTest', 'testNan', 138, 144]

[875, 'src/test/java', 'org.apache.commons.math3.linear', 'RealVectorFormatAbstractTest', 'testPositiveInfinity', 146, 154]

[875, 'src/test/java', 'org.apache.commons.math3.linear', 'RealVectorFormatAbstractTest', 'tesNegativeInfinity', 156, 164]

    @Test
    public void testNan() {
        ArrayRealVector c = new ArrayRealVector(new double[] {Double.NaN, Double.NaN, Double.NaN});
        String expected = "{(NaN); (NaN); (NaN)}";
        String actual = realVectorFormat.format(c);
        Assert.assertEquals(expected, actual);
    }

    @Test
    public void testPositiveInfinity() {
        ArrayRealVector c = new ArrayRealVector(new double[] {
                Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY
        });
        String expected = "{(Infinity); (Infinity); (Infinity)}";
        String actual = realVectorFormat.format(c);
        Assert.assertEquals(expected, actual);
    }

    @Test
    public void tesNegativeInfinity() {
        ArrayRealVector c = new ArrayRealVector(new double[] {
                Double.NEGATIVE_INFINITY, Double.NEGATIVE_INFINITY, Double.NEGATIVE_INFINITY
        });
        String expected = "{(-Infinity); (-Infinity); (-Infinity)}";
        String actual = realVectorFormat.format(c);
        Assert.assertEquals(expected, actual);
    }

[904, 'src/test/java', 'org.apache.commons.math3.linear', 'HessenbergTransformerTest', 'testNonSquare', 51, 59]

[904, 'src/test/java', 'org.apache.commons.math3.linear', 'LUDecompositionTest', 'testNonSquare', 71, 79]

[904, 'src/test/java', 'org.apache.commons.math3.linear', 'SchurTransformerTest', 'testNonSquare', 52, 60]

[904, 'src/test/java', 'org.apache.commons.math3.linear', 'TriDiagonalTransformerTest', 'testNonSquare', 43, 51]

    @Test
    public void testNonSquare() {
        try {
            new HessenbergTransformer(MatrixUtils.createRealMatrix(new double[3][2]));
            Assert.fail("an exception should have been thrown");
        } catch (NonSquareMatrixException ime) {
            // expected behavior
        }
    }

    @Test
    public void testNonSquare() {
        try {
            new LUDecomposition(MatrixUtils.createRealMatrix(new double[3][2]));
            Assert.fail("Expecting NonSquareMatrixException");
        } catch (NonSquareMatrixException ime) {
            // expected behavior
        }
    }

    @Test
    public void testNonSquare() {
        try {
            new SchurTransformer(MatrixUtils.createRealMatrix(new double[3][2]));
            Assert.fail("an exception should have been thrown");
        } catch (NonSquareMatrixException ime) {
            // expected behavior
        }
    }

    @Test
    public void testNonSquare() {
        try {
            new TriDiagonalTransformer(MatrixUtils.createRealMatrix(new double[3][2]));
            Assert.fail("an exception should have been thrown");
        } catch (NonSquareMatrixException ime) {
            // expected behavior
        }
    }

[931, 'src/test/java', 'org.apache.commons.math3.stat.descriptive.moment', 'MeanTest', 'testSmallSamples', 52, 58]

[931, 'src/test/java', 'org.apache.commons.math3.stat.descriptive.moment', 'StandardDeviationTest', 'testNaN', 53, 59]

[931, 'src/test/java', 'org.apache.commons.math3.stat.descriptive.moment', 'VarianceTest', 'testNaN', 58, 64]

    @Test
    public void testSmallSamples() {
        Mean mean = new Mean();
        Assert.assertTrue(Double.isNaN(mean.getResult()));
        mean.increment(1d);
        Assert.assertEquals(1d, mean.getResult(), 0);
    }

    @Test
    public void testNaN() {
        StandardDeviation std = new StandardDeviation();
        Assert.assertTrue(Double.isNaN(std.getResult()));
        std.increment(1d);
        Assert.assertEquals(0d, std.getResult(), 0);
    }

    @Test
    public void testNaN() {
        StandardDeviation std = new StandardDeviation();
        Assert.assertTrue(Double.isNaN(std.getResult()));
        std.increment(1d);
        Assert.assertEquals(0d, std.getResult(), 0);
    }

[986, 'src/test/java', 'org.apache.commons.math3.complex', 'ComplexTest', 'testScalarAdd', 131, 137]

[986, 'src/test/java', 'org.apache.commons.math3.complex', 'ComplexTest', 'testScalarDivide', 263, 269]

[986, 'src/test/java', 'org.apache.commons.math3.complex', 'ComplexTest', 'testScalarSubtract', 465, 471]

    @Test
    public void testScalarAdd() {
        Complex x = new Complex(3.0, 4.0);
        double yDouble = 2.0;
        Complex yComplex = new Complex(yDouble);
        Assert.assertEquals(x.add(yComplex), x.add(yDouble));
    }

    @Test
    public void testScalarDivide() {
        Complex x = new Complex(3.0, 4.0);
        double yDouble = 2.0;
        Complex yComplex = new Complex(yDouble);
        Assert.assertEquals(x.divide(yComplex), x.divide(yDouble));
    }

    @Test
    public void testScalarSubtract() {
        Complex x = new Complex(3.0, 4.0);
        double yDouble = 2.0;
        Complex yComplex = new Complex(yDouble);
        Assert.assertEquals(x.subtract(yComplex), x.subtract(yDouble));
    }

[1006, 'src/test/java', 'org.apache.commons.math3.analysis.solvers', 'MullerSolver2Test', 'testExpm1Function', 97, 120]

[1006, 'src/test/java', 'org.apache.commons.math3.analysis.solvers', 'MullerSolverTest', 'testExpm1Function', 99, 122]

[1006, 'src/test/java', 'org.apache.commons.math3.analysis.solvers', 'RiddersSolverTest', 'testExpm1Function', 93, 116]

    @Test
    public void testExpm1Function() {
        UnivariateFunction f = new Expm1();
        UnivariateSolver solver = new MullerSolver2();
        double min, max, expected, result, tolerance;

        min = -1.0; max = 2.0; expected = 0.0;
        tolerance = FastMath.max(solver.getAbsoluteAccuracy(),
                    FastMath.abs(expected * solver.getRelativeAccuracy()));
        result = solver.solve(100, f, min, max);
        Assert.assertEquals(expected, result, tolerance);

        min = -20.0; max = 10.0; expected = 0.0;
        tolerance = FastMath.max(solver.getAbsoluteAccuracy(),
                    FastMath.abs(expected * solver.getRelativeAccuracy()));
        result = solver.solve(100, f, min, max);
        Assert.assertEquals(expected, result, tolerance);

        min = -50.0; max = 100.0; expected = 0.0;
        tolerance = FastMath.max(solver.getAbsoluteAccuracy(),
                    FastMath.abs(expected * solver.getRelativeAccuracy()));
        result = solver.solve(100, f, min, max);
        Assert.assertEquals(expected, result, tolerance);
    }

    @Test
    public void testExpm1Function() {
        UnivariateFunction f = new Expm1();
        UnivariateSolver solver = new MullerSolver();
        double min, max, expected, result, tolerance;

        min = -1.0; max = 2.0; expected = 0.0;
        tolerance = FastMath.max(solver.getAbsoluteAccuracy(),
                    FastMath.abs(expected * solver.getRelativeAccuracy()));
        result = solver.solve(100, f, min, max);
        Assert.assertEquals(expected, result, tolerance);

        min = -20.0; max = 10.0; expected = 0.0;
        tolerance = FastMath.max(solver.getAbsoluteAccuracy(),
                    FastMath.abs(expected * solver.getRelativeAccuracy()));
        result = solver.solve(100, f, min, max);
        Assert.assertEquals(expected, result, tolerance);

        min = -50.0; max = 100.0; expected = 0.0;
        tolerance = FastMath.max(solver.getAbsoluteAccuracy(),
                    FastMath.abs(expected * solver.getRelativeAccuracy()));
        result = solver.solve(100, f, min, max);
        Assert.assertEquals(expected, result, tolerance);
    }

    @Test
    public void testExpm1Function() {
        UnivariateFunction f = new Expm1();
        UnivariateSolver solver = new RiddersSolver();
        double min, max, expected, result, tolerance;

        min = -1.0; max = 2.0; expected = 0.0;
        tolerance = FastMath.max(solver.getAbsoluteAccuracy(),
                    FastMath.abs(expected * solver.getRelativeAccuracy()));
        result = solver.solve(100, f, min, max);
        Assert.assertEquals(expected, result, tolerance);

        min = -20.0; max = 10.0; expected = 0.0;
        tolerance = FastMath.max(solver.getAbsoluteAccuracy(),
                    FastMath.abs(expected * solver.getRelativeAccuracy()));
        result = solver.solve(100, f, min, max);
        Assert.assertEquals(expected, result, tolerance);

        min = -50.0; max = 100.0; expected = 0.0;
        tolerance = FastMath.max(solver.getAbsoluteAccuracy(),
                    FastMath.abs(expected * solver.getRelativeAccuracy()));
        result = solver.solve(100, f, min, max);
        Assert.assertEquals(expected, result, tolerance);
    }

[1030, 'src/test/java', 'org.apache.commons.math3.linear', 'QRDecompositionTest', 'testAEqualQR', 81, 98]

[1030, 'src/test/java', 'org.apache.commons.math3.linear', 'QRDecompositionTest', 'testQOrthogonal', 107, 124]

[1030, 'src/test/java', 'org.apache.commons.math3.linear', 'RRQRDecompositionTest', 'testAPEqualQR', 80, 97]

[1030, 'src/test/java', 'org.apache.commons.math3.linear', 'RRQRDecompositionTest', 'testQOrthogonal', 106, 123]

    @Test
    public void testAEqualQR() {
        checkAEqualQR(MatrixUtils.createRealMatrix(testData3x3NonSingular));

        checkAEqualQR(MatrixUtils.createRealMatrix(testData3x3Singular));

        checkAEqualQR(MatrixUtils.createRealMatrix(testData3x4));

        checkAEqualQR(MatrixUtils.createRealMatrix(testData4x3));

        Random r = new Random(643895747384642l);
        int    p = (5 * BlockRealMatrix.BLOCK_SIZE) / 4;
        int    q = (7 * BlockRealMatrix.BLOCK_SIZE) / 4;
        checkAEqualQR(createTestMatrix(r, p, q));

        checkAEqualQR(createTestMatrix(r, q, p));

    }

    @Test
    public void testQOrthogonal() {
        checkQOrthogonal(MatrixUtils.createRealMatrix(testData3x3NonSingular));

        checkQOrthogonal(MatrixUtils.createRealMatrix(testData3x3Singular));

        checkQOrthogonal(MatrixUtils.createRealMatrix(testData3x4));

        checkQOrthogonal(MatrixUtils.createRealMatrix(testData4x3));

        Random r = new Random(643895747384642l);
        int    p = (5 * BlockRealMatrix.BLOCK_SIZE) / 4;
        int    q = (7 * BlockRealMatrix.BLOCK_SIZE) / 4;
        checkQOrthogonal(createTestMatrix(r, p, q));

        checkQOrthogonal(createTestMatrix(r, q, p));

    }

    @Test
    public void testAPEqualQR() {
        checkAPEqualQR(MatrixUtils.createRealMatrix(testData3x3NonSingular));

        checkAPEqualQR(MatrixUtils.createRealMatrix(testData3x3Singular));

        checkAPEqualQR(MatrixUtils.createRealMatrix(testData3x4));

        checkAPEqualQR(MatrixUtils.createRealMatrix(testData4x3));

        Random r = new Random(643895747384642l);
        int    p = (5 * BlockRealMatrix.BLOCK_SIZE) / 4;
        int    q = (7 * BlockRealMatrix.BLOCK_SIZE) / 4;
        checkAPEqualQR(createTestMatrix(r, p, q));

        checkAPEqualQR(createTestMatrix(r, q, p));

    }

    @Test
    public void testQOrthogonal() {
        checkQOrthogonal(MatrixUtils.createRealMatrix(testData3x3NonSingular));

        checkQOrthogonal(MatrixUtils.createRealMatrix(testData3x3Singular));

        checkQOrthogonal(MatrixUtils.createRealMatrix(testData3x4));

        checkQOrthogonal(MatrixUtils.createRealMatrix(testData4x3));

        Random r = new Random(643895747384642l);
        int    p = (5 * BlockRealMatrix.BLOCK_SIZE) / 4;
        int    q = (7 * BlockRealMatrix.BLOCK_SIZE) / 4;
        checkQOrthogonal(createTestMatrix(r, p, q));

        checkQOrthogonal(createTestMatrix(r, q, p));

    }

[1091, 'src/test/java', 'org.apache.commons.math3.linear', 'RealVectorFormatAbstractTest', 'testSimpleWithDecimals', 53, 63]

[1091, 'src/test/java', 'org.apache.commons.math3.linear', 'RealVectorFormatAbstractTest', 'testSimpleWithDecimalsTrunc', 65, 75]

[1091, 'src/test/java', 'org.apache.commons.math3.linear', 'RealVectorFormatAbstractTest', 'testParseSimpleWithDecimals', 187, 197]

[1091, 'src/test/java', 'org.apache.commons.math3.linear', 'RealVectorFormatAbstractTest', 'testParseSimpleWithDecimalsTrunc', 199, 209]

[1091, 'src/test/java', 'org.apache.commons.math3.linear', 'RealVectorFormatAbstractTest', 'testParseNonDefaultSetting', 271, 281]

    @Test
    public void testSimpleWithDecimals() {
        ArrayRealVector c = new ArrayRealVector(new double[] {1.23, 1.43, 1.63});
        String expected =
            "{1"    + getDecimalCharacter() +
            "23; 1" + getDecimalCharacter() +
            "43; 1" + getDecimalCharacter() +
            "63}";
        String actual = realVectorFormat.format(c);
        Assert.assertEquals(expected, actual);
    }

    @Test
    public void testSimpleWithDecimalsTrunc() {
        ArrayRealVector c = new ArrayRealVector(new double[] {1.232323232323, 1.43434343434343, 1.633333333333});
        String expected =
            "{1"    + getDecimalCharacter() +
            "2323232323; 1" + getDecimalCharacter() +
            "4343434343; 1" + getDecimalCharacter() +
            "6333333333}";
        String actual = realVectorFormat.format(c);
        Assert.assertEquals(expected, actual);
    }

    @Test
    public void testParseSimpleWithDecimals() {
        String source =
            "{1" + getDecimalCharacter() +
            "23; 1" + getDecimalCharacter() +
            "43; 1" + getDecimalCharacter() +
            "63}";
        ArrayRealVector expected = new ArrayRealVector(new double[] {1.23, 1.43, 1.63});
        ArrayRealVector actual = realVectorFormat.parse(source);
        Assert.assertEquals(expected, actual);
    }

    @Test
    public void testParseSimpleWithDecimalsTrunc() {
        String source =
            "{1" + getDecimalCharacter() +
            "2323; 1" + getDecimalCharacter() +
            "4343; 1" + getDecimalCharacter() +
            "6333}";
        ArrayRealVector expected = new ArrayRealVector(new double[] {1.2323, 1.4343, 1.6333});
        ArrayRealVector actual = realVectorFormat.parse(source);
        Assert.assertEquals(expected, actual);
    }

    @Test
    public void testParseNonDefaultSetting() {
        String source =
            "[1" + getDecimalCharacter() +
            "2323 : 1" + getDecimalCharacter() +
            "4343 : 1" + getDecimalCharacter() +
            "6333]";
        ArrayRealVector expected = new ArrayRealVector(new double[] {1.2323, 1.4343, 1.6333});
        ArrayRealVector actual = realVectorFormatSquare.parse(source);
        Assert.assertEquals(expected, actual);
    }

[1113, 'src/test/java', 'org.apache.commons.math3.analysis.integration', 'RombergIntegratorTest', 'testSinFunction', 42, 61]

[1113, 'src/test/java', 'org.apache.commons.math3.analysis.integration', 'SimpsonIntegratorTest', 'testSinFunction', 41, 60]

[1113, 'src/test/java', 'org.apache.commons.math3.analysis.integration', 'TrapezoidIntegratorTest', 'testSinFunction', 41, 60]

    @Test
    public void testSinFunction() {
        UnivariateFunction f = new Sin();
        UnivariateIntegrator integrator = new RombergIntegrator();
        double min, max, expected, result, tolerance;

        min = 0; max = FastMath.PI; expected = 2;
        tolerance = FastMath.abs(expected * integrator.getRelativeAccuracy());
        result = integrator.integrate(100, f, min, max);
        Assert.assertTrue(integrator.getEvaluations() < 50);
        Assert.assertTrue(integrator.getIterations()  < 10);
        Assert.assertEquals(expected, result, tolerance);

        min = -FastMath.PI/3; max = 0; expected = -0.5;
        tolerance = FastMath.abs(expected * integrator.getRelativeAccuracy());
        result = integrator.integrate(100, f, min, max);
        Assert.assertTrue(integrator.getEvaluations() < 50);
        Assert.assertTrue(integrator.getIterations()  < 10);
        Assert.assertEquals(expected, result, tolerance);
    }

    @Test
    public void testSinFunction() {
        UnivariateFunction f = new Sin();
        UnivariateIntegrator integrator = new SimpsonIntegrator();
        double min, max, expected, result, tolerance;

        min = 0; max = FastMath.PI; expected = 2;
        tolerance = FastMath.abs(expected * integrator.getRelativeAccuracy());
        result = integrator.integrate(1000, f, min, max);
        Assert.assertTrue(integrator.getEvaluations() < 100);
        Assert.assertTrue(integrator.getIterations()  < 10);
        Assert.assertEquals(expected, result, tolerance);

        min = -FastMath.PI/3; max = 0; expected = -0.5;
        tolerance = FastMath.abs(expected * integrator.getRelativeAccuracy());
        result = integrator.integrate(1000, f, min, max);
        Assert.assertTrue(integrator.getEvaluations() < 50);
        Assert.assertTrue(integrator.getIterations()  < 10);
        Assert.assertEquals(expected, result, tolerance);
    }

    @Test
    public void testSinFunction() {
        UnivariateFunction f = new Sin();
        UnivariateIntegrator integrator = new TrapezoidIntegrator();
        double min, max, expected, result, tolerance;

        min = 0; max = FastMath.PI; expected = 2;
        tolerance = FastMath.abs(expected * integrator.getRelativeAccuracy());
        result = integrator.integrate(10000, f, min, max);
        Assert.assertTrue(integrator.getEvaluations() < 2500);
        Assert.assertTrue(integrator.getIterations()  < 15);
        Assert.assertEquals(expected, result, tolerance);

        min = -FastMath.PI/3; max = 0; expected = -0.5;
        tolerance = FastMath.abs(expected * integrator.getRelativeAccuracy());
        result = integrator.integrate(10000, f, min, max);
        Assert.assertTrue(integrator.getEvaluations() < 2500);
        Assert.assertTrue(integrator.getIterations()  < 15);
        Assert.assertEquals(expected, result, tolerance);
    }

[1168, 'src/test/java', 'org.apache.commons.math3.linear', 'BlockFieldMatrixTest', 'testSetRowMatrix', 762, 781]

[1168, 'src/test/java', 'org.apache.commons.math3.linear', 'BlockFieldMatrixTest', 'testSetColumnMatrix', 826, 845]

[1168, 'src/test/java', 'org.apache.commons.math3.linear', 'FieldMatrixImplTest', 'testSetRowMatrix', 604, 623]

[1168, 'src/test/java', 'org.apache.commons.math3.linear', 'FieldMatrixImplTest', 'testSetColumnMatrix', 648, 667]

    @Test
    public void testSetRowMatrix() {
        FieldMatrix<Fraction> m = new BlockFieldMatrix<Fraction>(subTestData);
        FieldMatrix<Fraction> mRow3 = new BlockFieldMatrix<Fraction>(subRow3);
        Assert.assertNotSame(mRow3, m.getRowMatrix(0));
        m.setRowMatrix(0, mRow3);
        Assert.assertEquals(mRow3, m.getRowMatrix(0));
        try {
            m.setRowMatrix(-1, mRow3);
            Assert.fail("Expecting OutOfRangeException");
        } catch (OutOfRangeException ex) {
            // expected
        }
        try {
            m.setRowMatrix(0, m);
            Assert.fail("Expecting MatrixDimensionMismatchException");
        } catch (MatrixDimensionMismatchException ex) {
            // expected
        }
    }

    @Test
    public void testSetColumnMatrix() {
        FieldMatrix<Fraction> m = new BlockFieldMatrix<Fraction>(subTestData);
        FieldMatrix<Fraction> mColumn3 = new BlockFieldMatrix<Fraction>(subColumn3);
        Assert.assertNotSame(mColumn3, m.getColumnMatrix(1));
        m.setColumnMatrix(1, mColumn3);
        Assert.assertEquals(mColumn3, m.getColumnMatrix(1));
        try {
            m.setColumnMatrix(-1, mColumn3);
            Assert.fail("Expecting OutOfRangeException");
        } catch (OutOfRangeException ex) {
            // expected
        }
        try {
            m.setColumnMatrix(0, m);
            Assert.fail("Expecting MatrixDimensionMismatchException");
        } catch (MatrixDimensionMismatchException ex) {
            // expected
        }
    }

    @Test
    public void testSetRowMatrix() {
        FieldMatrix<Fraction> m = new Array2DRowFieldMatrix<Fraction>(subTestData);
        FieldMatrix<Fraction> mRow3 = new Array2DRowFieldMatrix<Fraction>(subRow3);
        Assert.assertNotSame(mRow3, m.getRowMatrix(0));
        m.setRowMatrix(0, mRow3);
        Assert.assertEquals(mRow3, m.getRowMatrix(0));
        try {
            m.setRowMatrix(-1, mRow3);
            Assert.fail("Expecting OutOfRangeException");
        } catch (OutOfRangeException ex) {
            // expected
        }
        try {
            m.setRowMatrix(0, m);
            Assert.fail("Expecting MatrixDimensionMismatchException");
        } catch (MatrixDimensionMismatchException ex) {
            // expected
        }
    }

    @Test
    public void testSetColumnMatrix() {
        FieldMatrix<Fraction> m = new Array2DRowFieldMatrix<Fraction>(subTestData);
        FieldMatrix<Fraction> mColumn3 = new Array2DRowFieldMatrix<Fraction>(subColumn3);
        Assert.assertNotSame(mColumn3, m.getColumnMatrix(1));
        m.setColumnMatrix(1, mColumn3);
        Assert.assertEquals(mColumn3, m.getColumnMatrix(1));
        try {
            m.setColumnMatrix(-1, mColumn3);
            Assert.fail("Expecting OutOfRangeException");
        } catch (OutOfRangeException ex) {
            // expected
        }
        try {
            m.setColumnMatrix(0, m);
            Assert.fail("Expecting MatrixDimensionMismatchException");
        } catch (MatrixDimensionMismatchException ex) {
            // expected
        }
    }

[1249, 'src/test/java', 'org.apache.commons.math3.geometry.euclidean.threed', 'SubLineTest', 'testIntersectionInsideOutside', 124, 130]

[1249, 'src/test/java', 'org.apache.commons.math3.geometry.euclidean.threed', 'SubLineTest', 'testIntersectionBoundaryOutside', 140, 146]

[1249, 'src/test/java', 'org.apache.commons.math3.geometry.euclidean.threed', 'SubLineTest', 'testIntersectionOutsideOutside', 148, 154]

[1249, 'src/test/java', 'org.apache.commons.math3.geometry.euclidean.threed', 'SubLineTest', 'testIntersectionNotIntersecting', 156, 162]

    @Test
    public void testIntersectionInsideOutside() throws MathIllegalArgumentException {
        SubLine sub1 = new SubLine(new Vector3D(1, 1, 1), new Vector3D(3, 1, 1), 1.0e-10);
        SubLine sub2 = new SubLine(new Vector3D(2, 0, 0), new Vector3D(2, 0.5, 0.5), 1.0e-10);
        Assert.assertNull(sub1.intersection(sub2, true));
        Assert.assertNull(sub1.intersection(sub2, false));
    }

    @Test
    public void testIntersectionBoundaryOutside() throws MathIllegalArgumentException {
        SubLine sub1 = new SubLine(new Vector3D(1, 1, 1), new Vector3D(2, 1, 1), 1.0e-10);
        SubLine sub2 = new SubLine(new Vector3D(2, 0, 0), new Vector3D(2, 0.5, 0.5), 1.0e-10);
        Assert.assertNull(sub1.intersection(sub2, true));
        Assert.assertNull(sub1.intersection(sub2, false));
    }

    @Test
    public void testIntersectionOutsideOutside() throws MathIllegalArgumentException {
        SubLine sub1 = new SubLine(new Vector3D(1, 1, 1), new Vector3D(1.5, 1, 1), 1.0e-10);
        SubLine sub2 = new SubLine(new Vector3D(2, 0, 0), new Vector3D(2, 0.5, 0.5), 1.0e-10);
        Assert.assertNull(sub1.intersection(sub2, true));
        Assert.assertNull(sub1.intersection(sub2, false));
    }

    @Test
    public void testIntersectionNotIntersecting() throws MathIllegalArgumentException {
        SubLine sub1 = new SubLine(new Vector3D(1, 1, 1), new Vector3D(1.5, 1, 1), 1.0e-10);
        SubLine sub2 = new SubLine(new Vector3D(2, 3, 0), new Vector3D(2, 3, 0.5), 1.0e-10);
        Assert.assertNull(sub1.intersection(sub2, true));
        Assert.assertNull(sub1.intersection(sub2, false));
    }

[1269, 'src/test/java', 'org.apache.commons.math3.distribution', 'ChiSquaredDistributionTest', 'testDfAccessors', 97, 101]

[1269, 'src/test/java', 'org.apache.commons.math3.distribution', 'ExponentialDistributionTest', 'testMeanAccessors', 111, 115]

[1269, 'src/test/java', 'org.apache.commons.math3.distribution', 'TDistributionTest', 'testDfAccessors', 116, 120]

    @Test
    public void testDfAccessors() {
        ChiSquaredDistribution distribution = (ChiSquaredDistribution) getDistribution();
        Assert.assertEquals(5d, distribution.getDegreesOfFreedom(), Double.MIN_VALUE);
    }

    @Test
    public void testMeanAccessors() {
        ExponentialDistribution distribution = (ExponentialDistribution) getDistribution();
        Assert.assertEquals(5d, distribution.getMean(), Double.MIN_VALUE);
    }

    @Test
    public void testDfAccessors() {
        TDistribution dist = (TDistribution) getDistribution();
        Assert.assertEquals(5d, dist.getDegreesOfFreedom(), Double.MIN_VALUE);
    }

[1282, 'src/test/java', 'org.apache.commons.math3.optim.nonlinear.scalar.noderiv', 'SimplexOptimizerMultiDirectionalTest', 'testMinimize2', 69, 88]

[1282, 'src/test/java', 'org.apache.commons.math3.optim.nonlinear.scalar.noderiv', 'SimplexOptimizerNelderMeadTest', 'testMinimize2', 74, 93]

[1282, 'src/test/java', 'org.apache.commons.math3.optim.nonlinear.scalar.noderiv', 'SimplexOptimizerNelderMeadTest', 'testMaximize2', 116, 135]

    @Test
    public void testMinimize2() {
        SimplexOptimizer optimizer = new SimplexOptimizer(1e-11, 1e-30);
        final FourExtrema fourExtrema = new FourExtrema();

        final PointValuePair optimum
            = optimizer.optimize(new MaxEval(200),
                                 new ObjectiveFunction(fourExtrema),
                                 GoalType.MINIMIZE,
                                 new InitialGuess(new double[] { 1, 0 }),
                                 new MultiDirectionalSimplex(new double[] { 0.2, 0.2 }));
        Assert.assertEquals(fourExtrema.xP, optimum.getPoint()[0], 2e-8);
        Assert.assertEquals(fourExtrema.yM, optimum.getPoint()[1], 3e-6);
        Assert.assertEquals(fourExtrema.valueXpYm, optimum.getValue(), 2e-12);
        Assert.assertTrue(optimizer.getEvaluations() > 120);
        Assert.assertTrue(optimizer.getEvaluations() < 150);

        // Check that the number of iterations is updated (MATH-949).
        Assert.assertTrue(optimizer.getIterations() > 0);
    }

    @Test
    public void testMinimize2() {
        SimplexOptimizer optimizer = new SimplexOptimizer(1e-10, 1e-30);
        final FourExtrema fourExtrema = new FourExtrema();

        final PointValuePair optimum
            = optimizer.optimize(new MaxEval(100),
                                 new ObjectiveFunction(fourExtrema),
                                 GoalType.MINIMIZE,
                                 new InitialGuess(new double[] { 1, 0 }),
                                 new NelderMeadSimplex(new double[] { 0.2, 0.2 }));
        Assert.assertEquals(fourExtrema.xP, optimum.getPoint()[0], 5e-6);
        Assert.assertEquals(fourExtrema.yM, optimum.getPoint()[1], 6e-6);
        Assert.assertEquals(fourExtrema.valueXpYm, optimum.getValue(), 1e-11);
        Assert.assertTrue(optimizer.getEvaluations() > 60);
        Assert.assertTrue(optimizer.getEvaluations() < 90);

        // Check that the number of iterations is updated (MATH-949).
        Assert.assertTrue(optimizer.getIterations() > 0);
    }

    @Test
    public void testMaximize2() {
        SimplexOptimizer optimizer = new SimplexOptimizer(1e-10, 1e-30);
        final FourExtrema fourExtrema = new FourExtrema();

        final PointValuePair optimum
            = optimizer.optimize(new MaxEval(100),
                                 new ObjectiveFunction(fourExtrema),
                                 GoalType.MAXIMIZE,
                                 new InitialGuess(new double[] { 1, 0 }),
                                 new NelderMeadSimplex(new double[] { 0.2, 0.2 }));
        Assert.assertEquals(fourExtrema.xP, optimum.getPoint()[0], 4e-6);
        Assert.assertEquals(fourExtrema.yP, optimum.getPoint()[1], 5e-6);
        Assert.assertEquals(fourExtrema.valueXpYp, optimum.getValue(), 7e-12);
        Assert.assertTrue(optimizer.getEvaluations() > 60);
        Assert.assertTrue(optimizer.getEvaluations() < 90);

        // Check that the number of iterations is updated (MATH-949).
        Assert.assertTrue(optimizer.getIterations() > 0);
    }

[1319, 'src/test/java', 'org.apache.commons.math3.linear', 'QRDecompositionTest', 'testRUpperTriangular', 134, 157]

[1319, 'src/test/java', 'org.apache.commons.math3.linear', 'QRDecompositionTest', 'testHTrapezoidal', 171, 194]

[1319, 'src/test/java', 'org.apache.commons.math3.linear', 'RRQRDecompositionTest', 'testRUpperTriangular', 133, 156]

[1319, 'src/test/java', 'org.apache.commons.math3.linear', 'RRQRDecompositionTest', 'testHTrapezoidal', 170, 193]

    @Test
    public void testRUpperTriangular() {
        RealMatrix matrix = MatrixUtils.createRealMatrix(testData3x3NonSingular);
        checkUpperTriangular(new QRDecomposition(matrix).getR());

        matrix = MatrixUtils.createRealMatrix(testData3x3Singular);
        checkUpperTriangular(new QRDecomposition(matrix).getR());

        matrix = MatrixUtils.createRealMatrix(testData3x4);
        checkUpperTriangular(new QRDecomposition(matrix).getR());

        matrix = MatrixUtils.createRealMatrix(testData4x3);
        checkUpperTriangular(new QRDecomposition(matrix).getR());

        Random r = new Random(643895747384642l);
        int    p = (5 * BlockRealMatrix.BLOCK_SIZE) / 4;
        int    q = (7 * BlockRealMatrix.BLOCK_SIZE) / 4;
        matrix = createTestMatrix(r, p, q);
        checkUpperTriangular(new QRDecomposition(matrix).getR());

        matrix = createTestMatrix(r, p, q);
        checkUpperTriangular(new QRDecomposition(matrix).getR());

    }

    @Test
    public void testHTrapezoidal() {
        RealMatrix matrix = MatrixUtils.createRealMatrix(testData3x3NonSingular);
        checkTrapezoidal(new QRDecomposition(matrix).getH());

        matrix = MatrixUtils.createRealMatrix(testData3x3Singular);
        checkTrapezoidal(new QRDecomposition(matrix).getH());

        matrix = MatrixUtils.createRealMatrix(testData3x4);
        checkTrapezoidal(new QRDecomposition(matrix).getH());

        matrix = MatrixUtils.createRealMatrix(testData4x3);
        checkTrapezoidal(new QRDecomposition(matrix).getH());

        Random r = new Random(643895747384642l);
        int    p = (5 * BlockRealMatrix.BLOCK_SIZE) / 4;
        int    q = (7 * BlockRealMatrix.BLOCK_SIZE) / 4;
        matrix = createTestMatrix(r, p, q);
        checkTrapezoidal(new QRDecomposition(matrix).getH());

        matrix = createTestMatrix(r, p, q);
        checkTrapezoidal(new QRDecomposition(matrix).getH());

    }

    @Test
    public void testRUpperTriangular() {
        RealMatrix matrix = MatrixUtils.createRealMatrix(testData3x3NonSingular);
        checkUpperTriangular(new RRQRDecomposition(matrix).getR());

        matrix = MatrixUtils.createRealMatrix(testData3x3Singular);
        checkUpperTriangular(new RRQRDecomposition(matrix).getR());

        matrix = MatrixUtils.createRealMatrix(testData3x4);
        checkUpperTriangular(new RRQRDecomposition(matrix).getR());

        matrix = MatrixUtils.createRealMatrix(testData4x3);
        checkUpperTriangular(new RRQRDecomposition(matrix).getR());

        Random r = new Random(643895747384642l);
        int    p = (5 * BlockRealMatrix.BLOCK_SIZE) / 4;
        int    q = (7 * BlockRealMatrix.BLOCK_SIZE) / 4;
        matrix = createTestMatrix(r, p, q);
        checkUpperTriangular(new RRQRDecomposition(matrix).getR());

        matrix = createTestMatrix(r, p, q);
        checkUpperTriangular(new RRQRDecomposition(matrix).getR());

    }

    @Test
    public void testHTrapezoidal() {
        RealMatrix matrix = MatrixUtils.createRealMatrix(testData3x3NonSingular);
        checkTrapezoidal(new RRQRDecomposition(matrix).getH());

        matrix = MatrixUtils.createRealMatrix(testData3x3Singular);
        checkTrapezoidal(new RRQRDecomposition(matrix).getH());

        matrix = MatrixUtils.createRealMatrix(testData3x4);
        checkTrapezoidal(new RRQRDecomposition(matrix).getH());

        matrix = MatrixUtils.createRealMatrix(testData4x3);
        checkTrapezoidal(new RRQRDecomposition(matrix).getH());

        Random r = new Random(643895747384642l);
        int    p = (5 * BlockRealMatrix.BLOCK_SIZE) / 4;
        int    q = (7 * BlockRealMatrix.BLOCK_SIZE) / 4;
        matrix = createTestMatrix(r, p, q);
        checkTrapezoidal(new RRQRDecomposition(matrix).getH());

        matrix = createTestMatrix(r, p, q);
        checkTrapezoidal(new RRQRDecomposition(matrix).getH());

    }

[1325, 'src/test/java', 'org.apache.commons.math3.distribution', 'LogNormalDistributionTest', 'testDensity', 182, 193]

[1325, 'src/test/java', 'org.apache.commons.math3.distribution', 'NormalDistributionTest', 'testDensity', 140, 147]

[1325, 'src/test/java', 'org.apache.commons.math3.distribution', 'ParetoDistributionTest', 'testDensity', 153, 160]

    @Test
    public void testDensity() {
        double [] x = new double[]{-2, -1, 0, 1, 2};
        // R 2.13: print(dlnorm(c(-2,-1,0,1,2)), digits=10)
        checkDensity(0, 1, x, new double[] { 0.0000000000, 0.0000000000,
                                             0.0000000000, 0.3989422804,
                                             0.1568740193 });
        // R 2.13: print(dlnorm(c(-2,-1,0,1,2), mean=1.1), digits=10)
        checkDensity(1.1, 1, x, new double[] { 0.0000000000, 0.0000000000,
                                               0.0000000000, 0.2178521770,
                                               0.1836267118});
    }

    @Test
    public void testDensity() {
        double [] x = new double[]{-2, -1, 0, 1, 2};
        // R 2.5: print(dnorm(c(-2,-1,0,1,2)), digits=10)
        checkDensity(0, 1, x, new double[]{0.05399096651, 0.24197072452, 0.39894228040, 0.24197072452, 0.05399096651});
        // R 2.5: print(dnorm(c(-2,-1,0,1,2), mean=1.1), digits=10)
        checkDensity(1.1, 1, x, new double[]{0.003266819056,0.043983595980,0.217852177033,0.396952547477,0.266085249899});
    }

    @Test
    public void testDensity() {
        double [] x = new double[]{-2, -1, 0, 1, 2};
        // R 2.14: print(dpareto(c(-2,-1,0,1,2), scale=1, shape=1), digits=10)
        checkDensity(1, 1, x, new double[] { 0.00, 0.00, 0.00, 1.00, 0.25 });
        // R 2.14: print(dpareto(c(-2,-1,0,1,2), scale=1.1, shape=1), digits=10)
        checkDensity(1.1, 1, x, new double[] { 0.000, 0.000, 0.000, 0.000, 0.275 });
    }

[1341, 'src/test/java', 'org.apache.commons.math3.fitting', 'GaussianCurveFitterTest', 'testFit04', 255, 263]

[1341, 'src/test/java', 'org.apache.commons.math3.fitting', 'GaussianCurveFitterTest', 'testFit05', 268, 276]

[1341, 'src/test/java', 'org.apache.commons.math3.fitting', 'GaussianCurveFitterTest', 'testFit06', 281, 289]

    @Test
    public void testFit04() {
        GaussianCurveFitter fitter = GaussianCurveFitter.create();
        double[] parameters = fitter.fit(createDataset(DATASET2).toList());

        Assert.assertEquals(233003.2967252038, parameters[0], 1e-4);
        Assert.assertEquals(-10.654887521095983, parameters[1], 1e-4);
        Assert.assertEquals(4.335937353196641, parameters[2], 1e-4);
    }

    @Test
    public void testFit05() {
        GaussianCurveFitter fitter = GaussianCurveFitter.create();
        double[] parameters = fitter.fit(createDataset(DATASET3).toList());

        Assert.assertEquals(283863.81929180305, parameters[0], 1e-4);
        Assert.assertEquals(-13.29641995105174, parameters[1], 1e-4);
        Assert.assertEquals(1.7297330293549908, parameters[2], 1e-4);
    }

    @Test
    public void testFit06() {
        GaussianCurveFitter fitter = GaussianCurveFitter.create();
        double[] parameters = fitter.fit(createDataset(DATASET4).toList());

        Assert.assertEquals(285250.66754309234, parameters[0], 1e-4);
        Assert.assertEquals(-13.528375695228455, parameters[1], 1e-4);
        Assert.assertEquals(1.5204344894331614, parameters[2], 1e-4);
    }

[1345, 'src/test/java', 'org.apache.commons.math3.geometry.euclidean.twod.hull', 'ConvexHullGenerator2DAbstractTest', 'testCollinearPoints', 117, 128]

[1345, 'src/test/java', 'org.apache.commons.math3.geometry.euclidean.twod.hull', 'ConvexHullGenerator2DAbstractTest', 'testCollinearPointsReverse', 130, 141]

[1345, 'src/test/java', 'org.apache.commons.math3.geometry.euclidean.twod.hull', 'ConvexHullGenerator2DAbstractTest', 'testIdenticalPoints', 169, 180]

    @Test
    public void testCollinearPoints() {
        final Collection<Vector2D> points = new ArrayList<Vector2D>();
        points.add(new Vector2D(1, 1));
        points.add(new Vector2D(2, 2));
        points.add(new Vector2D(2, 4));
        points.add(new Vector2D(4, 1));
        points.add(new Vector2D(10, 1));

        final ConvexHull2D hull = generator.generate(points);
        checkConvexHull(points, hull);
    }

    @Test
    public void testCollinearPointsReverse() {
        final Collection<Vector2D> points = new ArrayList<Vector2D>();
        points.add(new Vector2D(1, 1));
        points.add(new Vector2D(2, 2));
        points.add(new Vector2D(2, 4));
        points.add(new Vector2D(10, 1));
        points.add(new Vector2D(4, 1));

        final ConvexHull2D hull = generator.generate(points);
        checkConvexHull(points, hull);
    }

    @Test
    public void testIdenticalPoints() {
        final Collection<Vector2D> points = new ArrayList<Vector2D>();
        points.add(new Vector2D(1, 1));
        points.add(new Vector2D(2, 2));
        points.add(new Vector2D(2, 4));
        points.add(new Vector2D(4, 1));
        points.add(new Vector2D(1, 1));

        final ConvexHull2D hull = generator.generate(points);
        checkConvexHull(points, hull);
    }

[1401, 'src/test/java', 'org.apache.commons.math3.linear', 'DiagonalMatrixTest', 'testSetNonDiagonalZero', 281, 286]

[1401, 'src/test/java', 'org.apache.commons.math3.linear', 'DiagonalMatrixTest', 'testAddNonDiagonalZero', 294, 299]

[1401, 'src/test/java', 'org.apache.commons.math3.linear', 'DiagonalMatrixTest', 'testMultiplyNonDiagonalEntry', 301, 306]

[1401, 'src/test/java', 'org.apache.commons.math3.linear', 'DiagonalMatrixTest', 'testMultiplyNonDiagonalZero', 308, 313]

    @Test
    public void testSetNonDiagonalZero() {
        final DiagonalMatrix diag = new DiagonalMatrix(3);
        diag.setEntry(1, 2, 0.0);
        Assert.assertEquals(0.0, diag.getEntry(1, 2), Precision.SAFE_MIN);
    }

    @Test
    public void testAddNonDiagonalZero() {
        final DiagonalMatrix diag = new DiagonalMatrix(3);
        diag.addToEntry(1, 2, 0.0);
        Assert.assertEquals(0.0, diag.getEntry(1, 2), Precision.SAFE_MIN);
    }

    @Test
    public void testMultiplyNonDiagonalEntry() {
        final DiagonalMatrix diag = new DiagonalMatrix(3);
        diag.multiplyEntry(1, 2, 3.4);
        Assert.assertEquals(0.0, diag.getEntry(1, 2), Precision.SAFE_MIN);
    }

    @Test
    public void testMultiplyNonDiagonalZero() {
        final DiagonalMatrix diag = new DiagonalMatrix(3);
        diag.multiplyEntry(1, 2, 0.0);
        Assert.assertEquals(0.0, diag.getEntry(1, 2), Precision.SAFE_MIN);
    }

[1407, 'src/test/java', 'org.apache.commons.math3.distribution', 'GumbelDistributionTest', 'testSupport', 35, 41]

[1407, 'src/test/java', 'org.apache.commons.math3.distribution', 'LaplaceDistributionTest', 'testSupport', 35, 41]

[1407, 'src/test/java', 'org.apache.commons.math3.distribution', 'LogisticsDistributionTest', 'testSupport', 35, 41]

    @Test
    public void testSupport() {
        GumbelDistribution d = makeDistribution();
        Assert.assertTrue(Double.isInfinite(d.getSupportLowerBound()));
        Assert.assertTrue(Double.isInfinite(d.getSupportUpperBound()));
        Assert.assertTrue(d.isSupportConnected());
    }

    @Test
    public void testSupport() {
        LaplaceDistribution d = makeDistribution();
        Assert.assertTrue(Double.isInfinite(d.getSupportLowerBound()));
        Assert.assertTrue(Double.isInfinite(d.getSupportUpperBound()));
        Assert.assertTrue(d.isSupportConnected());
    }

    @Test
    public void testSupport() {
        LogisticDistribution d = makeDistribution();
        Assert.assertTrue(Double.isInfinite(d.getSupportLowerBound()));
        Assert.assertTrue(Double.isInfinite(d.getSupportUpperBound()));
        Assert.assertTrue(d.isSupportConnected());
    }

[1432, 'src/test/java', 'org.apache.commons.math3.special', 'ErfTest', 'testErf1960', 37, 49]

[1432, 'src/test/java', 'org.apache.commons.math3.special', 'ErfTest', 'testErf2807', 65, 77]

[1432, 'src/test/java', 'org.apache.commons.math3.special', 'ErfTest', 'testErf3291', 79, 91]

    @Test
    public void testErf1960() {
        double x = 1.960 / FastMath.sqrt(2.0);
        double actual = Erf.erf(x);
        double expected = 0.95;
        Assert.assertEquals(expected, actual, 1.0e-5);
        Assert.assertEquals(1 - actual, Erf.erfc(x), 1.0e-15);

        actual = Erf.erf(-x);
        expected = -expected;
        Assert.assertEquals(expected, actual, 1.0e-5);
        Assert.assertEquals(1 - actual, Erf.erfc(-x), 1.0e-15);
    }

    @Test
    public void testErf2807() {
        double x = 2.807 / FastMath.sqrt(2.0);
        double actual = Erf.erf(x);
        double expected = 0.995;
        Assert.assertEquals(expected, actual, 1.0e-5);
        Assert.assertEquals(1 - actual, Erf.erfc(x), 1.0e-15);

        actual = Erf.erf(-x);
        expected = -expected;
        Assert.assertEquals(expected, actual, 1.0e-5);
        Assert.assertEquals(1 - actual, Erf.erfc(-x), 1.0e-15);
    }

    @Test
    public void testErf3291() {
        double x = 3.291 / FastMath.sqrt(2.0);
        double actual = Erf.erf(x);
        double expected = 0.999;
        Assert.assertEquals(expected, actual, 1.0e-5);
        Assert.assertEquals(1 - expected, Erf.erfc(x), 1.0e-5);

        actual = Erf.erf(-x);
        expected = -expected;
        Assert.assertEquals(expected, actual, 1.0e-5);
        Assert.assertEquals(1 - expected, Erf.erfc(-x), 1.0e-5);
    }

[1470, 'src/test/java', 'org.apache.commons.math3.linear', 'RealVectorFormatAbstractTest', 'testParseNan', 283, 288]

[1470, 'src/test/java', 'org.apache.commons.math3.linear', 'RealVectorFormatAbstractTest', 'testParsePositiveInfinity', 290, 297]

[1470, 'src/test/java', 'org.apache.commons.math3.linear', 'RealVectorFormatAbstractTest', 'testParseNegativeInfinity', 299, 306]

    @Test
    public void testParseNan() {
        String source = "{(NaN); (NaN); (NaN)}";
        ArrayRealVector actual = realVectorFormat.parse(source);
        Assert.assertEquals(new ArrayRealVector(new double[] {Double.NaN, Double.NaN, Double.NaN}), actual);
    }

    @Test
    public void testParsePositiveInfinity() {
        String source = "{(Infinity); (Infinity); (Infinity)}";
        ArrayRealVector actual = realVectorFormat.parse(source);
        Assert.assertEquals(new ArrayRealVector(new double[] {
                Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY
        }), actual);
    }

    @Test
    public void testParseNegativeInfinity() {
        String source = "{(-Infinity); (-Infinity); (-Infinity)}";
        ArrayRealVector actual = realVectorFormat.parse(source);
        Assert.assertEquals(new ArrayRealVector(new double[] {
                Double.NEGATIVE_INFINITY, Double.NEGATIVE_INFINITY, Double.NEGATIVE_INFINITY
        }), actual);
    }

[1471, 'src/test/java', 'org.apache.commons.math3.fraction', 'FractionFormatTest', 'testParseOne1', 139, 146]

[1471, 'src/test/java', 'org.apache.commons.math3.fraction', 'FractionFormatTest', 'testParseOne2', 148, 155]

[1471, 'src/test/java', 'org.apache.commons.math3.fraction', 'FractionFormatTest', 'testParseZero1', 157, 164]

    @Test
    public void testParseOne1() {
        String source = "1 / 1";
        Fraction c = properFormat.parse(source);
        Assert.assertNotNull(c);
        Assert.assertEquals(1, c.getNumerator());
        Assert.assertEquals(1, c.getDenominator());
    }

    @Test
    public void testParseOne2() {
        String source = "10 / 10";
        Fraction c = properFormat.parse(source);
        Assert.assertNotNull(c);
        Assert.assertEquals(1, c.getNumerator());
        Assert.assertEquals(1, c.getDenominator());
    }

    @Test
    public void testParseZero1() {
        String source = "0 / 1";
        Fraction c = properFormat.parse(source);
        Assert.assertNotNull(c);
        Assert.assertEquals(0, c.getNumerator());
        Assert.assertEquals(1, c.getDenominator());
    }

[1490, 'src/test/java', 'org.apache.commons.math3.linear', 'QRSolverTest', 'testSolveDimensionErrors', 72, 89]

[1490, 'src/test/java', 'org.apache.commons.math3.linear', 'QRSolverTest', 'testSolveRankErrors', 92, 109]

[1490, 'src/test/java', 'org.apache.commons.math3.linear', 'RRQRSolverTest', 'testSolveDimensionErrors', 72, 90]

    @Test
    public void testSolveDimensionErrors() {
        DecompositionSolver solver =
            new QRDecomposition(MatrixUtils.createRealMatrix(testData3x3NonSingular)).getSolver();
        RealMatrix b = MatrixUtils.createRealMatrix(new double[2][2]);
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
    public void testSolveRankErrors() {
        DecompositionSolver solver =
            new QRDecomposition(MatrixUtils.createRealMatrix(testData3x3Singular)).getSolver();
        RealMatrix b = MatrixUtils.createRealMatrix(new double[3][2]);
        try {
            solver.solve(b);
            Assert.fail("an exception should have been thrown");
        } catch (SingularMatrixException iae) {
            // expected behavior
        }
        try {
            solver.solve(b.getColumnVector(0));
            Assert.fail("an exception should have been thrown");
        } catch (SingularMatrixException iae) {
            // expected behavior
        }
    }

    @Test
    public void testSolveDimensionErrors() {
        DecompositionSolver solver =
            new RRQRDecomposition(MatrixUtils.createRealMatrix(testData3x3NonSingular)).getSolver();
        RealMatrix b = MatrixUtils.createRealMatrix(new double[2][2]);
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

[1503, 'src/test/java', 'org.apache.commons.math3.linear', 'HessenbergTransformerTest', 'testPOrthogonal', 68, 72]

[1503, 'src/test/java', 'org.apache.commons.math3.linear', 'HessenbergTransformerTest', 'testPTOrthogonal', 74, 78]

[1503, 'src/test/java', 'org.apache.commons.math3.linear', 'HessenbergTransformerTest', 'testHessenbergForm', 80, 84]

[1503, 'src/test/java', 'org.apache.commons.math3.linear', 'TriDiagonalTransformerTest', 'testQOrthogonal', 89, 93]

[1503, 'src/test/java', 'org.apache.commons.math3.linear', 'TriDiagonalTransformerTest', 'testQTOrthogonal', 95, 99]

[1503, 'src/test/java', 'org.apache.commons.math3.linear', 'TriDiagonalTransformerTest', 'testTTriDiagonal', 107, 111]

    @Test
    public void testPOrthogonal() {
        checkOrthogonal(new HessenbergTransformer(MatrixUtils.createRealMatrix(testSquare5)).getP());
        checkOrthogonal(new HessenbergTransformer(MatrixUtils.createRealMatrix(testSquare3)).getP());
    }

    @Test
    public void testPTOrthogonal() {
        checkOrthogonal(new HessenbergTransformer(MatrixUtils.createRealMatrix(testSquare5)).getPT());
        checkOrthogonal(new HessenbergTransformer(MatrixUtils.createRealMatrix(testSquare3)).getPT());
    }

    @Test
    public void testHessenbergForm() {
        checkHessenbergForm(new HessenbergTransformer(MatrixUtils.createRealMatrix(testSquare5)).getH());
        checkHessenbergForm(new HessenbergTransformer(MatrixUtils.createRealMatrix(testSquare3)).getH());
    }

    @Test
    public void testQOrthogonal() {
        checkOrthogonal(new TriDiagonalTransformer(MatrixUtils.createRealMatrix(testSquare5)).getQ());
        checkOrthogonal(new TriDiagonalTransformer(MatrixUtils.createRealMatrix(testSquare3)).getQ());
    }

    @Test
    public void testQTOrthogonal() {
        checkOrthogonal(new TriDiagonalTransformer(MatrixUtils.createRealMatrix(testSquare5)).getQT());
        checkOrthogonal(new TriDiagonalTransformer(MatrixUtils.createRealMatrix(testSquare3)).getQT());
    }

    @Test
    public void testTTriDiagonal() {
        checkTriDiagonal(new TriDiagonalTransformer(MatrixUtils.createRealMatrix(testSquare5)).getT());
        checkTriDiagonal(new TriDiagonalTransformer(MatrixUtils.createRealMatrix(testSquare3)).getT());
    }

[1607, 'src/test/java', 'org.apache.commons.math3.ode.nonstiff', 'ClassicalRungeKuttaStepInterpolatorTest', 'derivativesConsistency', 41, 49]

[1607, 'src/test/java', 'org.apache.commons.math3.ode.nonstiff', 'EulerStepInterpolatorTest', 'derivativesConsistency', 132, 140]

[1607, 'src/test/java', 'org.apache.commons.math3.ode.nonstiff', 'GillStepInterpolatorTest', 'testDerivativesConsistency', 41, 49]

[1607, 'src/test/java', 'org.apache.commons.math3.ode.nonstiff', 'LutherStepInterpolatorTest', 'derivativesConsistency', 41, 49]

[1607, 'src/test/java', 'org.apache.commons.math3.ode.nonstiff', 'MidpointStepInterpolatorTest', 'testDerivativesConsistency', 42, 50]

[1607, 'src/test/java', 'org.apache.commons.math3.ode.nonstiff', 'ThreeEighthesStepInterpolatorTest', 'derivativesConsistency', 41, 49]

  @Test
  public void derivativesConsistency()
      throws DimensionMismatchException, NumberIsTooSmallException,
             MaxCountExceededException, NoBracketingException {
    TestProblem3 pb = new TestProblem3();
    double step = (pb.getFinalTime() - pb.getInitialTime()) * 0.001;
    ClassicalRungeKuttaIntegrator integ = new ClassicalRungeKuttaIntegrator(step);
    StepInterpolatorTestUtils.checkDerivativesConsistency(integ, pb, 0.01, 6.6e-12);
  }

  @Test
  public void derivativesConsistency()
      throws DimensionMismatchException, NumberIsTooSmallException,
             MaxCountExceededException, NoBracketingException {
    TestProblem3 pb = new TestProblem3();
    double step = (pb.getFinalTime() - pb.getInitialTime()) * 0.001;
    EulerIntegrator integ = new EulerIntegrator(step);
    StepInterpolatorTestUtils.checkDerivativesConsistency(integ, pb, 0.01, 5.1e-12);
  }

  @Test
  public void testDerivativesConsistency()
      throws DimensionMismatchException, NumberIsTooSmallException,
             MaxCountExceededException, NoBracketingException {
    TestProblem3 pb = new TestProblem3();
    double step = (pb.getFinalTime() - pb.getInitialTime()) * 0.001;
    GillIntegrator integ = new GillIntegrator(step);
    StepInterpolatorTestUtils.checkDerivativesConsistency(integ, pb, 0.01, 6.6e-12);
  }

    @Test
    public void derivativesConsistency()
            throws DimensionMismatchException, NumberIsTooSmallException,
            MaxCountExceededException, NoBracketingException {
        TestProblem3 pb = new TestProblem3();
        double step = (pb.getFinalTime() - pb.getInitialTime()) * 0.001;
        LutherIntegrator integ = new LutherIntegrator(step);
        StepInterpolatorTestUtils.checkDerivativesConsistency(integ, pb, 0.01, 6.5e-12);
    }

  @Test
  public void testDerivativesConsistency()
      throws DimensionMismatchException, NumberIsTooSmallException,
             MaxCountExceededException, NoBracketingException {
    TestProblem3 pb = new TestProblem3();
    double step = (pb.getFinalTime() - pb.getInitialTime()) * 0.001;
    MidpointIntegrator integ = new MidpointIntegrator(step);
    StepInterpolatorTestUtils.checkDerivativesConsistency(integ, pb, 0.01, 6.6e-12);
  }

  @Test
  public void derivativesConsistency()
      throws DimensionMismatchException, NumberIsTooSmallException,
             MaxCountExceededException, NoBracketingException {
    TestProblem3 pb = new TestProblem3();
    double step = (pb.getFinalTime() - pb.getInitialTime()) * 0.001;
    ThreeEighthesIntegrator integ = new ThreeEighthesIntegrator(step);
    StepInterpolatorTestUtils.checkDerivativesConsistency(integ, pb, 0.01, 6.6e-12);
  }

[1639, 'src/test/java', 'org.apache.commons.math3.stat.regression', 'OLSMultipleLinearRegressionTest', 'testNoDataNPECalculateBeta', 803, 807]

[1639, 'src/test/java', 'org.apache.commons.math3.stat.regression', 'OLSMultipleLinearRegressionTest', 'testNoDataNPECalculateHat', 809, 813]

[1639, 'src/test/java', 'org.apache.commons.math3.stat.regression', 'OLSMultipleLinearRegressionTest', 'testNoDataNPESSTO', 815, 819]

    @Test(expected=NullPointerException.class)
    public void testNoDataNPECalculateBeta() {
        OLSMultipleLinearRegression model = new OLSMultipleLinearRegression();
        model.calculateBeta();
    }

    @Test(expected=NullPointerException.class)
    public void testNoDataNPECalculateHat() {
        OLSMultipleLinearRegression model = new OLSMultipleLinearRegression();
        model.calculateHat();
    }

    @Test(expected=NullPointerException.class)
    public void testNoDataNPESSTO() {
        OLSMultipleLinearRegression model = new OLSMultipleLinearRegression();
        model.calculateTotalSumOfSquares();
    }

[1677, 'src/test/java', 'org.apache.commons.math3.ode.nonstiff', 'ClassicalRungeKuttaIntegratorTest', 'testKepler', 244, 257]

[1677, 'src/test/java', 'org.apache.commons.math3.ode.nonstiff', 'GillIntegratorTest', 'testKepler', 164, 177]

[1677, 'src/test/java', 'org.apache.commons.math3.ode.nonstiff', 'LutherIntegratorTest', 'testKepler', 244, 257]

[1677, 'src/test/java', 'org.apache.commons.math3.ode.nonstiff', 'ThreeEighthesIntegratorTest', 'testKepler', 164, 177]

  @Test
  public void testKepler()
      throws DimensionMismatchException, NumberIsTooSmallException,
             MaxCountExceededException, NoBracketingException {

    final TestProblem3 pb  = new TestProblem3(0.9);
    double step = (pb.getFinalTime() - pb.getInitialTime()) * 0.0003;

    FirstOrderIntegrator integ = new ClassicalRungeKuttaIntegrator(step);
    integ.addStepHandler(new KeplerHandler(pb));
    integ.integrate(pb,
                    pb.getInitialTime(), pb.getInitialState(),
                    pb.getFinalTime(), new double[pb.getDimension()]);
  }

  @Test
  public void testKepler()
      throws DimensionMismatchException, NumberIsTooSmallException,
             MaxCountExceededException, NoBracketingException {

    final TestProblem3 pb  = new TestProblem3(0.9);
    double step = (pb.getFinalTime() - pb.getInitialTime()) * 0.0003;

    FirstOrderIntegrator integ = new GillIntegrator(step);
    integ.addStepHandler(new KeplerStepHandler(pb));
    integ.integrate(pb,
                    pb.getInitialTime(), pb.getInitialState(),
                    pb.getFinalTime(), new double[pb.getDimension()]);
  }

    @Test
    public void testKepler()
            throws DimensionMismatchException, NumberIsTooSmallException,
            MaxCountExceededException, NoBracketingException {

        final TestProblem3 pb  = new TestProblem3(0.9);
        double step = (pb.getFinalTime() - pb.getInitialTime()) * 0.0003;

        FirstOrderIntegrator integ = new LutherIntegrator(step);
        integ.addStepHandler(new KeplerHandler(pb));
        integ.integrate(pb,
                        pb.getInitialTime(), pb.getInitialState(),
                        pb.getFinalTime(), new double[pb.getDimension()]);
    }

  @Test
  public void testKepler()
      throws DimensionMismatchException, NumberIsTooSmallException,
             MaxCountExceededException, NoBracketingException {

    final TestProblem3 pb  = new TestProblem3(0.9);
    double step = (pb.getFinalTime() - pb.getInitialTime()) * 0.0003;

    FirstOrderIntegrator integ = new ThreeEighthesIntegrator(step);
    integ.addStepHandler(new KeplerHandler(pb));
    integ.integrate(pb,
                    pb.getInitialTime(), pb.getInitialState(),
                    pb.getFinalTime(), new double[pb.getDimension()]);
  }

[1688, 'src/test/java', 'org.apache.commons.math3.geometry.euclidean.twod.hull', 'ConvexHullGenerator2DAbstractTest', 'testCollinearPointsIncluded', 143, 154]

[1688, 'src/test/java', 'org.apache.commons.math3.geometry.euclidean.twod.hull', 'ConvexHullGenerator2DAbstractTest', 'testCollinearPointsIncludedReverse', 156, 167]

[1688, 'src/test/java', 'org.apache.commons.math3.geometry.euclidean.twod.hull', 'ConvexHullGenerator2DAbstractTest', 'testIdenticalPoints2', 182, 193]

    @Test
    public void testCollinearPointsIncluded() {
        final Collection<Vector2D> points = new ArrayList<Vector2D>();
        points.add(new Vector2D(1, 1));
        points.add(new Vector2D(2, 2));
        points.add(new Vector2D(2, 4));
        points.add(new Vector2D(4, 1));
        points.add(new Vector2D(10, 1));

        final ConvexHull2D hull = createConvexHullGenerator(true).generate(points);
        checkConvexHull(points, hull, true);
    }

    @Test
    public void testCollinearPointsIncludedReverse() {
        final Collection<Vector2D> points = new ArrayList<Vector2D>();
        points.add(new Vector2D(1, 1));
        points.add(new Vector2D(2, 2));
        points.add(new Vector2D(2, 4));
        points.add(new Vector2D(10, 1));
        points.add(new Vector2D(4, 1));

        final ConvexHull2D hull = createConvexHullGenerator(true).generate(points);
        checkConvexHull(points, hull, true);
    }

    @Test
    public void testIdenticalPoints2() {
        final Collection<Vector2D> points = new ArrayList<Vector2D>();
        points.add(new Vector2D(1, 1));
        points.add(new Vector2D(2, 2));
        points.add(new Vector2D(2, 4));
        points.add(new Vector2D(4, 1));
        points.add(new Vector2D(1, 1));

        final ConvexHull2D hull = createConvexHullGenerator(true).generate(points);
        checkConvexHull(points, hull, true);
    }

[1710, 'src/test/java', 'org.apache.commons.math3.geometry.euclidean.threed', 'FieldVector3DTest', 'testDistance', 282, 302]

[1710, 'src/test/java', 'org.apache.commons.math3.geometry.euclidean.threed', 'FieldVector3DTest', 'testDistanceSq', 304, 324]

[1710, 'src/test/java', 'org.apache.commons.math3.geometry.euclidean.threed', 'FieldVector3DTest', 'testDistance1', 260, 280]

    @Test
    public void testDistance() {
        FieldVector3D<DerivativeStructure> v1 = createVector(1, -2, 3, 3);
        FieldVector3D<DerivativeStructure> v2 = createVector(-4, 2, 0, 3);
        Assert.assertEquals(0.0, FieldVector3D.distance(createVector(-1, 0, 0, 3), createVector(-1, 0, 0, 3)).getReal(), 0);
        DerivativeStructure distance = FieldVector3D.distance(v1, v2);
        Assert.assertEquals(FastMath.sqrt(50), distance.getReal(), 1.0e-12);
        Assert.assertEquals(0, distance.getPartialDerivative(1, 0, 0), 1.0e-12);
        Assert.assertEquals(0, distance.getPartialDerivative(0, 1, 0), 1.0e-12);
        Assert.assertEquals(0, distance.getPartialDerivative(0, 0, 1), 1.0e-12);
        distance = FieldVector3D.distance(v1, new Vector3D(-4, 2, 0));
        Assert.assertEquals(FastMath.sqrt(50), distance.getReal(), 1.0e-12);
        Assert.assertEquals( 5 / FastMath.sqrt(50), distance.getPartialDerivative(1, 0, 0), 1.0e-12);
        Assert.assertEquals(-4 / FastMath.sqrt(50), distance.getPartialDerivative(0, 1, 0), 1.0e-12);
        Assert.assertEquals( 3 / FastMath.sqrt(50), distance.getPartialDerivative(0, 0, 1), 1.0e-12);
        distance = FieldVector3D.distance(new Vector3D(-4, 2, 0), v1);
        Assert.assertEquals(FastMath.sqrt(50), distance.getReal(), 1.0e-12);
        Assert.assertEquals( 5 / FastMath.sqrt(50), distance.getPartialDerivative(1, 0, 0), 1.0e-12);
        Assert.assertEquals(-4 / FastMath.sqrt(50), distance.getPartialDerivative(0, 1, 0), 1.0e-12);
        Assert.assertEquals( 3 / FastMath.sqrt(50), distance.getPartialDerivative(0, 0, 1), 1.0e-12);
    }

    @Test
    public void testDistanceSq() {
        FieldVector3D<DerivativeStructure> v1 = createVector(1, -2, 3, 3);
        FieldVector3D<DerivativeStructure> v2 = createVector(-4, 2, 0, 3);
        Assert.assertEquals(0.0, FieldVector3D.distanceSq(createVector(-1, 0, 0, 3), createVector(-1, 0, 0, 3)).getReal(), 0);
        DerivativeStructure distanceSq = FieldVector3D.distanceSq(v1, v2);
        Assert.assertEquals(50.0, distanceSq.getReal(), 1.0e-12);
        Assert.assertEquals(0, distanceSq.getPartialDerivative(1, 0, 0), 1.0e-12);
        Assert.assertEquals(0, distanceSq.getPartialDerivative(0, 1, 0), 1.0e-12);
        Assert.assertEquals(0, distanceSq.getPartialDerivative(0, 0, 1), 1.0e-12);
        distanceSq = FieldVector3D.distanceSq(v1, new Vector3D(-4, 2, 0));
        Assert.assertEquals(50.0, distanceSq.getReal(), 1.0e-12);
        Assert.assertEquals(10, distanceSq.getPartialDerivative(1, 0, 0), 1.0e-12);
        Assert.assertEquals(-8, distanceSq.getPartialDerivative(0, 1, 0), 1.0e-12);
        Assert.assertEquals( 6, distanceSq.getPartialDerivative(0, 0, 1), 1.0e-12);
        distanceSq = FieldVector3D.distanceSq(new Vector3D(-4, 2, 0), v1);
        Assert.assertEquals(50.0, distanceSq.getReal(), 1.0e-12);
        Assert.assertEquals(10, distanceSq.getPartialDerivative(1, 0, 0), 1.0e-12);
        Assert.assertEquals(-8, distanceSq.getPartialDerivative(0, 1, 0), 1.0e-12);
        Assert.assertEquals( 6, distanceSq.getPartialDerivative(0, 0, 1), 1.0e-12);
  }

    @Test
    public void testDistance1() {
        FieldVector3D<DerivativeStructure> v1 = createVector(1, -2, 3, 3);
        FieldVector3D<DerivativeStructure> v2 = createVector(-4, 2, 0, 3);
        Assert.assertEquals(0.0, FieldVector3D.distance1(createVector(-1, 0, 0, 3), createVector(-1, 0, 0, 3)).getReal(), 0);
        DerivativeStructure distance = FieldVector3D.distance1(v1, v2);
        Assert.assertEquals(12.0, distance.getReal(), 1.0e-12);
        Assert.assertEquals(0, distance.getPartialDerivative(1, 0, 0), 1.0e-12);
        Assert.assertEquals(0, distance.getPartialDerivative(0, 1, 0), 1.0e-12);
        Assert.assertEquals(0, distance.getPartialDerivative(0, 0, 1), 1.0e-12);
        distance = FieldVector3D.distance1(v1, new Vector3D(-4, 2, 0));
        Assert.assertEquals(12.0, distance.getReal(), 1.0e-12);
        Assert.assertEquals( 1, distance.getPartialDerivative(1, 0, 0), 1.0e-12);
        Assert.assertEquals(-1, distance.getPartialDerivative(0, 1, 0), 1.0e-12);
        Assert.assertEquals( 1, distance.getPartialDerivative(0, 0, 1), 1.0e-12);
        distance = FieldVector3D.distance1(new Vector3D(-4, 2, 0), v1);
        Assert.assertEquals(12.0, distance.getReal(), 1.0e-12);
        Assert.assertEquals( 1, distance.getPartialDerivative(1, 0, 0), 1.0e-12);
        Assert.assertEquals(-1, distance.getPartialDerivative(0, 1, 0), 1.0e-12);
        Assert.assertEquals( 1, distance.getPartialDerivative(0, 0, 1), 1.0e-12);
    }

[1779, 'src/test/java', 'org.apache.commons.math3.linear', 'BiDiagonalTransformerTest', 'testDimensions', 41, 46]

[1779, 'src/test/java', 'org.apache.commons.math3.linear', 'BiDiagonalTransformerTest', 'testAEqualUSVt', 61, 66]

[1779, 'src/test/java', 'org.apache.commons.math3.linear', 'SingularValueDecompositionTest', 'testAEqualUSVt', 130, 135]

    @Test
    public void testDimensions() {
        checkdimensions(MatrixUtils.createRealMatrix(testSquare));
        checkdimensions(MatrixUtils.createRealMatrix(testNonSquare));
        checkdimensions(MatrixUtils.createRealMatrix(testNonSquare).transpose());
    }

    @Test
    public void testAEqualUSVt() {
        checkAEqualUSVt(MatrixUtils.createRealMatrix(testSquare));
        checkAEqualUSVt(MatrixUtils.createRealMatrix(testNonSquare));
        checkAEqualUSVt(MatrixUtils.createRealMatrix(testNonSquare).transpose());
    }

    @Test
    public void testAEqualUSVt() {
        checkAEqualUSVt(MatrixUtils.createRealMatrix(testSquare));
        checkAEqualUSVt(MatrixUtils.createRealMatrix(testNonSquare));
        checkAEqualUSVt(MatrixUtils.createRealMatrix(testNonSquare).transpose());
    }

[1797, 'src/test/java', 'org.apache.commons.math3.fitting.leastsquares', 'GaussNewtonOptimizerWithCholeskyTest', 'testMaxEvaluations', 77, 99]

[1797, 'src/test/java', 'org.apache.commons.math3.fitting.leastsquares', 'GaussNewtonOptimizerWithLUTest', 'testMaxEvaluations', 77, 99]

[1797, 'src/test/java', 'org.apache.commons.math3.fitting.leastsquares', 'GaussNewtonOptimizerWithQRTest', 'testMaxEvaluations', 63, 85]

[1797, 'src/test/java', 'org.apache.commons.math3.fitting.leastsquares', 'GaussNewtonOptimizerWithSVDTest', 'testMaxEvaluations', 54, 76]

    @Test
    public void testMaxEvaluations() throws Exception {
        try{
        CircleVectorial circle = new CircleVectorial();
        circle.addPoint( 30.0,  68.0);
        circle.addPoint( 50.0,  -6.0);
        circle.addPoint(110.0, -20.0);
        circle.addPoint( 35.0,  15.0);
        circle.addPoint( 45.0,  97.0);

        LeastSquaresProblem lsp = builder(circle)
                .checkerPair(new SimpleVectorValueChecker(1e-30, 1e-30))
                .maxIterations(Integer.MAX_VALUE)
                .start(new double[]{98.680, 47.345})
                .build();

        optimizer.optimize(lsp);

            fail(optimizer);
        }catch (TooManyEvaluationsException e){
            //expected
        }
    }

    @Test
    public void testMaxEvaluations() throws Exception {
        try{
        CircleVectorial circle = new CircleVectorial();
        circle.addPoint( 30.0,  68.0);
        circle.addPoint( 50.0,  -6.0);
        circle.addPoint(110.0, -20.0);
        circle.addPoint( 35.0,  15.0);
        circle.addPoint( 45.0,  97.0);

        LeastSquaresProblem lsp = builder(circle)
                .checkerPair(new SimpleVectorValueChecker(1e-30, 1e-30))
                .maxIterations(Integer.MAX_VALUE)
                .start(new double[]{98.680, 47.345})
                .build();

        optimizer.optimize(lsp);

            fail(optimizer);
        }catch (TooManyEvaluationsException e){
            //expected
        }
    }

    @Test
    public void testMaxEvaluations() throws Exception {
        try{
        CircleVectorial circle = new CircleVectorial();
        circle.addPoint( 30.0,  68.0);
        circle.addPoint( 50.0,  -6.0);
        circle.addPoint(110.0, -20.0);
        circle.addPoint( 35.0,  15.0);
        circle.addPoint( 45.0,  97.0);

        LeastSquaresProblem lsp = builder(circle)
                .checkerPair(new SimpleVectorValueChecker(1e-30, 1e-30))
                .maxIterations(Integer.MAX_VALUE)
                .start(new double[]{98.680, 47.345})
                .build();

        optimizer.optimize(lsp);

            fail(optimizer);
        }catch (TooManyEvaluationsException e){
            //expected
        }
    }

    @Test
    public void testMaxEvaluations() throws Exception {
        try{
        CircleVectorial circle = new CircleVectorial();
        circle.addPoint( 30.0,  68.0);
        circle.addPoint( 50.0,  -6.0);
        circle.addPoint(110.0, -20.0);
        circle.addPoint( 35.0,  15.0);
        circle.addPoint( 45.0,  97.0);

        LeastSquaresProblem lsp = builder(circle)
                .checkerPair(new SimpleVectorValueChecker(1e-30, 1e-30))
                .maxIterations(Integer.MAX_VALUE)
                .start(new double[]{98.680, 47.345})
                .build();

        optimizer.optimize(lsp);

            fail(optimizer);
        }catch (TooManyEvaluationsException e){
            //expected
        }
    }

[1829, 'src/test/java', 'org.apache.commons.math3.linear', 'RealMatrixFormatAbstractTest', 'testNegativeComponent', 85, 99]

[1829, 'src/test/java', 'org.apache.commons.math3.linear', 'RealMatrixFormatAbstractTest', 'testNegativeComponent2', 101, 115]

[1829, 'src/test/java', 'org.apache.commons.math3.linear', 'RealMatrixFormatAbstractTest', 'testNegativeSecondRow', 117, 131]

    @Test
    public void testNegativeComponent() {
        RealMatrix m = MatrixUtils.createRealMatrix(new double[][] {{-1.232323232323, 1.43, 1.63},
                                                                    {2.46, 2.46, 2.66}});
        String expected =
                "{{-1"    + getDecimalCharacter() +
                "2323232323,1" + getDecimalCharacter() +
                "43,1" + getDecimalCharacter() +
                "63},{2" + getDecimalCharacter() +
                "46,2" + getDecimalCharacter() +
                "46,2" + getDecimalCharacter() +
                "66}}";
        String actual = realMatrixFormat.format(m);
        Assert.assertEquals(expected, actual);
    }

    @Test
    public void testNegativeComponent2() {
        RealMatrix m = MatrixUtils.createRealMatrix(new double[][] {{1.23, -1.434343434343, 1.63},
                                                                    {2.46, 2.46, 2.66}});
        String expected =
                "{{1"    + getDecimalCharacter() +
                "23,-1" + getDecimalCharacter() +
                "4343434343,1" + getDecimalCharacter() +
                "63},{2" + getDecimalCharacter() +
                "46,2" + getDecimalCharacter() +
                "46,2" + getDecimalCharacter() +
                "66}}";
        String actual = realMatrixFormat.format(m);
        Assert.assertEquals(expected, actual);
    }

    @Test
    public void testNegativeSecondRow() {
        RealMatrix m = MatrixUtils.createRealMatrix(new double[][] {{1.23, 1.43, 1.63},
                                                                    {-2.66666666666, 2.46, 2.66}});
        String expected =
                "{{1"    + getDecimalCharacter() +
                "23,1" + getDecimalCharacter() +
                "43,1" + getDecimalCharacter() +
                "63},{-2" + getDecimalCharacter() +
                "6666666667,2" + getDecimalCharacter() +
                "46,2" + getDecimalCharacter() +
                "66}}";
        String actual = realMatrixFormat.format(m);
        Assert.assertEquals(expected, actual);
    }

[1840, 'src/test/java', 'org.apache.commons.math3.ode.nonstiff', 'ClassicalRungeKuttaIntegratorTest', 'testStepSize', 289, 314]

[1840, 'src/test/java', 'org.apache.commons.math3.ode.nonstiff', 'EulerIntegratorTest', 'testStepSize', 166, 191]

[1840, 'src/test/java', 'org.apache.commons.math3.ode.nonstiff', 'GillIntegratorTest', 'testStepSize', 220, 245]

[1840, 'src/test/java', 'org.apache.commons.math3.ode.nonstiff', 'LutherIntegratorTest', 'testStepSize', 285, 310]

[1840, 'src/test/java', 'org.apache.commons.math3.ode.nonstiff', 'MidpointIntegratorTest', 'testStepSize', 166, 191]

[1840, 'src/test/java', 'org.apache.commons.math3.ode.nonstiff', 'ThreeEighthesIntegratorTest', 'testStepSize', 214, 239]

  @Test
  public void testStepSize()
      throws DimensionMismatchException, NumberIsTooSmallException,
             MaxCountExceededException, NoBracketingException {
      final double step = 1.23456;
      FirstOrderIntegrator integ = new ClassicalRungeKuttaIntegrator(step);
      integ.addStepHandler(new StepHandler() {
          public void handleStep(StepInterpolator interpolator, boolean isLast) {
              if (! isLast) {
                  Assert.assertEquals(step,
                               interpolator.getCurrentTime() - interpolator.getPreviousTime(),
                               1.0e-12);
              }
          }
          public void init(double t0, double[] y0, double t) {
          }
      });
      integ.integrate(new FirstOrderDifferentialEquations() {
          public void computeDerivatives(double t, double[] y, double[] dot) {
              dot[0] = 1.0;
          }
          public int getDimension() {
              return 1;
          }
      }, 0.0, new double[] { 0.0 }, 5.0, new double[1]);
  }

  @Test
  public void testStepSize()
      throws DimensionMismatchException, NumberIsTooSmallException,
             MaxCountExceededException, NoBracketingException {
      final double step = 1.23456;
      FirstOrderIntegrator integ = new EulerIntegrator(step);
      integ.addStepHandler(new StepHandler() {
        public void handleStep(StepInterpolator interpolator, boolean isLast) {
            if (! isLast) {
                Assert.assertEquals(step,
                             interpolator.getCurrentTime() - interpolator.getPreviousTime(),
                             1.0e-12);
            }
        }
        public void init(double t0, double[] y0, double t) {
        }
      });
      integ.integrate(new FirstOrderDifferentialEquations() {
                          public void computeDerivatives(double t, double[] y, double[] dot) {
                              dot[0] = 1.0;
                          }
                          public int getDimension() {
                              return 1;
                          }
                      }, 0.0, new double[] { 0.0 }, 5.0, new double[1]);
  }

  @Test
  public void testStepSize()
      throws DimensionMismatchException, NumberIsTooSmallException,
             MaxCountExceededException, NoBracketingException {
      final double step = 1.23456;
      FirstOrderIntegrator integ = new GillIntegrator(step);
      integ.addStepHandler(new StepHandler() {
          public void handleStep(StepInterpolator interpolator, boolean isLast) {
              if (! isLast) {
                  Assert.assertEquals(step,
                               interpolator.getCurrentTime() - interpolator.getPreviousTime(),
                               1.0e-12);
              }
          }
          public void init(double t0, double[] y0, double t) {
          }
      });
      integ.integrate(new FirstOrderDifferentialEquations() {
          public void computeDerivatives(double t, double[] y, double[] dot) {
              dot[0] = 1.0;
          }
          public int getDimension() {
              return 1;
          }
      }, 0.0, new double[] { 0.0 }, 5.0, new double[1]);
  }

    @Test
    public void testStepSize()
            throws DimensionMismatchException, NumberIsTooSmallException,
            MaxCountExceededException, NoBracketingException {
        final double step = 1.23456;
        FirstOrderIntegrator integ = new LutherIntegrator(step);
        integ.addStepHandler(new StepHandler() {
            public void handleStep(StepInterpolator interpolator, boolean isLast) {
                if (! isLast) {
                    Assert.assertEquals(step,
                                        interpolator.getCurrentTime() - interpolator.getPreviousTime(),
                                        1.0e-12);
                }
            }
            public void init(double t0, double[] y0, double t) {
            }
        });
        integ.integrate(new FirstOrderDifferentialEquations() {
            public void computeDerivatives(double t, double[] y, double[] dot) {
                dot[0] = 1.0;
            }
            public int getDimension() {
                return 1;
            }
        }, 0.0, new double[] { 0.0 }, 5.0, new double[1]);
    }

  @Test
  public void testStepSize()
      throws DimensionMismatchException, NumberIsTooSmallException,
             MaxCountExceededException, NoBracketingException {
      final double step = 1.23456;
      FirstOrderIntegrator integ = new MidpointIntegrator(step);
      integ.addStepHandler(new StepHandler() {
          public void handleStep(StepInterpolator interpolator, boolean isLast) {
              if (! isLast) {
                  Assert.assertEquals(step,
                               interpolator.getCurrentTime() - interpolator.getPreviousTime(),
                               1.0e-12);
              }
          }
          public void init(double t0, double[] y0, double t) {
          }
      });
      integ.integrate(new FirstOrderDifferentialEquations() {
          public void computeDerivatives(double t, double[] y, double[] dot) {
              dot[0] = 1.0;
          }
          public int getDimension() {
              return 1;
          }
      }, 0.0, new double[] { 0.0 }, 5.0, new double[1]);
  }

  @Test
  public void testStepSize()
      throws DimensionMismatchException, NumberIsTooSmallException,
             MaxCountExceededException, NoBracketingException {
      final double step = 1.23456;
      FirstOrderIntegrator integ = new ThreeEighthesIntegrator(step);
      integ.addStepHandler(new StepHandler() {
          public void handleStep(StepInterpolator interpolator, boolean isLast) {
              if (! isLast) {
                  Assert.assertEquals(step,
                               interpolator.getCurrentTime() - interpolator.getPreviousTime(),
                               1.0e-12);
              }
          }
          public void init(double t0, double[] y0, double t) {
          }
      });
      integ.integrate(new FirstOrderDifferentialEquations() {
          public void computeDerivatives(double t, double[] y, double[] dot) {
              dot[0] = 1.0;
          }
          public int getDimension() {
              return 1;
          }
      }, 0.0, new double[] { 0.0 }, 5.0, new double[1]);
  }

[1883, 'src/test/java', 'org.apache.commons.math3.fraction', 'BigFractionFormatTest', 'testFormat', 47, 57]

[1883, 'src/test/java', 'org.apache.commons.math3.fraction', 'BigFractionFormatTest', 'testFormatZero', 71, 81]

[1883, 'src/test/java', 'org.apache.commons.math3.fraction', 'FractionFormatTest', 'testFormat', 45, 55]

[1883, 'src/test/java', 'org.apache.commons.math3.fraction', 'FractionFormatTest', 'testFormatZero', 69, 79]

    @Test
    public void testFormat() {
        BigFraction c = new BigFraction(1, 2);
        String expected = "1 / 2";

        String actual = properFormat.format(c);
        Assert.assertEquals(expected, actual);

        actual = improperFormat.format(c);
        Assert.assertEquals(expected, actual);
    }

    @Test
    public void testFormatZero() {
        BigFraction c = new BigFraction(0, 1);
        String expected = "0 / 1";

        String actual = properFormat.format(c);
        Assert.assertEquals(expected, actual);

        actual = improperFormat.format(c);
        Assert.assertEquals(expected, actual);
    }

    @Test
    public void testFormat() {
        Fraction c = new Fraction(1, 2);
        String expected = "1 / 2";

        String actual = properFormat.format(c);
        Assert.assertEquals(expected, actual);

        actual = improperFormat.format(c);
        Assert.assertEquals(expected, actual);
    }

    @Test
    public void testFormatZero() {
        Fraction c = new Fraction(0, 1);
        String expected = "0 / 1";

        String actual = properFormat.format(c);
        Assert.assertEquals(expected, actual);

        actual = improperFormat.format(c);
        Assert.assertEquals(expected, actual);
    }

[1884, 'src/test/java', 'org.apache.commons.math3.geometry.euclidean.threed', 'Vector3DTest', 'testDistance1', 158, 165]

[1884, 'src/test/java', 'org.apache.commons.math3.geometry.euclidean.threed', 'Vector3DTest', 'testDistance', 167, 174]

[1884, 'src/test/java', 'org.apache.commons.math3.geometry.euclidean.threed', 'Vector3DTest', 'testDistanceSq', 176, 184]

[1884, 'src/test/java', 'org.apache.commons.math3.geometry.euclidean.threed', 'Vector3DTest', 'testDistanceInf', 186, 193]

    @Test
    public void testDistance1() {
        Vector3D v1 = new Vector3D(1, -2, 3);
        Vector3D v2 = new Vector3D(-4, 2, 0);
        Assert.assertEquals(0.0, Vector3D.distance1(Vector3D.MINUS_I, Vector3D.MINUS_I), 0);
        Assert.assertEquals(12.0, Vector3D.distance1(v1, v2), 1.0e-12);
        Assert.assertEquals(v1.subtract(v2).getNorm1(), Vector3D.distance1(v1, v2), 1.0e-12);
    }

    @Test
    public void testDistance() {
        Vector3D v1 = new Vector3D(1, -2, 3);
        Vector3D v2 = new Vector3D(-4, 2, 0);
        Assert.assertEquals(0.0, Vector3D.distance(Vector3D.MINUS_I, Vector3D.MINUS_I), 0);
        Assert.assertEquals(FastMath.sqrt(50), Vector3D.distance(v1, v2), 1.0e-12);
        Assert.assertEquals(v1.subtract(v2).getNorm(), Vector3D.distance(v1, v2), 1.0e-12);
    }

    @Test
    public void testDistanceSq() {
        Vector3D v1 = new Vector3D(1, -2, 3);
        Vector3D v2 = new Vector3D(-4, 2, 0);
        Assert.assertEquals(0.0, Vector3D.distanceSq(Vector3D.MINUS_I, Vector3D.MINUS_I), 0);
        Assert.assertEquals(50.0, Vector3D.distanceSq(v1, v2), 1.0e-12);
        Assert.assertEquals(Vector3D.distance(v1, v2) * Vector3D.distance(v1, v2),
                            Vector3D.distanceSq(v1, v2), 1.0e-12);
  }

    @Test
    public void testDistanceInf() {
        Vector3D v1 = new Vector3D(1, -2, 3);
        Vector3D v2 = new Vector3D(-4, 2, 0);
        Assert.assertEquals(0.0, Vector3D.distanceInf(Vector3D.MINUS_I, Vector3D.MINUS_I), 0);
        Assert.assertEquals(5.0, Vector3D.distanceInf(v1, v2), 1.0e-12);
        Assert.assertEquals(v1.subtract(v2).getNormInf(), Vector3D.distanceInf(v1, v2), 1.0e-12);
    }

[2027, 'src/test/java', 'org.apache.commons.math3.complex', 'ComplexTest', 'testCos', 736, 741]

[2027, 'src/test/java', 'org.apache.commons.math3.complex', 'ComplexTest', 'testCosh', 760, 765]

[2027, 'src/test/java', 'org.apache.commons.math3.complex', 'ComplexTest', 'testSinh', 1020, 1025]

    @Test
    public void testCos() {
        Complex z = new Complex(3, 4);
        Complex expected = new Complex(-27.03495, -3.851153);
        TestUtils.assertEquals(expected, z.cos(), 1.0e-5);
    }

    @Test
    public void testCosh() {
        Complex z = new Complex(3, 4);
        Complex expected = new Complex(-6.58066, -7.58155);
        TestUtils.assertEquals(expected, z.cosh(), 1.0e-5);
    }

    @Test
    public void testSinh() {
        Complex z = new Complex(3, 4);
        Complex expected = new Complex(-6.54812, -7.61923);
        TestUtils.assertEquals(expected, z.sinh(), 1.0e-5);
    }

[2079, 'src/test/java', 'org.apache.commons.math3.complex', 'ComplexTest', 'testAsin', 682, 687]

[2079, 'src/test/java', 'org.apache.commons.math3.complex', 'ComplexTest', 'testAtan', 707, 712]

[2079, 'src/test/java', 'org.apache.commons.math3.complex', 'ComplexTest', 'testLog', 842, 847]

    @Test
    public void testAsin() {
        Complex z = new Complex(3, 4);
        Complex expected = new Complex(0.633984, 2.30551);
        TestUtils.assertEquals(expected, z.asin(), 1.0e-5);
    }

    @Test
    public void testAtan() {
        Complex z = new Complex(3, 4);
        Complex expected = new Complex(1.44831, 0.158997);
        TestUtils.assertEquals(expected, z.atan(), 1.0e-5);
    }

    @Test
    public void testLog() {
        Complex z = new Complex(3, 4);
        Complex expected = new Complex(1.60944, 0.927295);
        TestUtils.assertEquals(expected, z.log(), 1.0e-5);
    }

[2123, 'src/test/java', 'org.apache.commons.math3.stat.interval', 'AgrestiCoullIntervalTest', 'testStandardInterval', 36, 41]

[2123, 'src/test/java', 'org.apache.commons.math3.stat.interval', 'ClopperPearsonIntervalTest', 'testStandardInterval', 36, 41]

[2123, 'src/test/java', 'org.apache.commons.math3.stat.interval', 'NormalApproximationIntervalTest', 'testStandardInterval', 36, 41]

[2123, 'src/test/java', 'org.apache.commons.math3.stat.interval', 'WilsonScoreIntervalTest', 'testStandardInterval', 36, 41]

    @Test
    public void testStandardInterval() {
        ConfidenceInterval confidenceInterval = createStandardTestInterval();
        Assert.assertEquals(0.07993521, confidenceInterval.getLowerBound(), 1E-5);
        Assert.assertEquals(0.1243704, confidenceInterval.getUpperBound(), 1E-5);
    }

    @Test
    public void testStandardInterval() {
        ConfidenceInterval confidenceInterval = createStandardTestInterval();
        Assert.assertEquals(0.07873857, confidenceInterval.getLowerBound(), 1E-5);
        Assert.assertEquals(0.1248658, confidenceInterval.getUpperBound(), 1E-5);
    }

    @Test
    public void testStandardInterval() {
        ConfidenceInterval confidenceInterval = createStandardTestInterval();
        Assert.assertEquals(0.07793197, confidenceInterval.getLowerBound(), 1E-5);
        Assert.assertEquals(0.1220680, confidenceInterval.getUpperBound(), 1E-5);
    }

    @Test
    public void testStandardInterval() {
        ConfidenceInterval confidenceInterval = createStandardTestInterval();
        Assert.assertEquals(0.08003919, confidenceInterval.getLowerBound(), 1E-5);
        Assert.assertEquals(0.1242664, confidenceInterval.getUpperBound(), 1E-5);
    }

[2137, 'src/test/java', 'org.apache.commons.math3.analysis.solvers', 'MullerSolver2Test', 'testQuinticFunction', 67, 90]

[2137, 'src/test/java', 'org.apache.commons.math3.analysis.solvers', 'MullerSolverTest', 'testQuinticFunction', 67, 90]

[2137, 'src/test/java', 'org.apache.commons.math3.analysis.solvers', 'RiddersSolverTest', 'testQuinticFunction', 65, 88]

    @Test
    public void testQuinticFunction() {
        UnivariateFunction f = new QuinticFunction();
        UnivariateSolver solver = new MullerSolver2();
        double min, max, expected, result, tolerance;

        min = -0.4; max = 0.2; expected = 0.0;
        tolerance = FastMath.max(solver.getAbsoluteAccuracy(),
                    FastMath.abs(expected * solver.getRelativeAccuracy()));
        result = solver.solve(100, f, min, max);
        Assert.assertEquals(expected, result, tolerance);

        min = 0.75; max = 1.5; expected = 1.0;
        tolerance = FastMath.max(solver.getAbsoluteAccuracy(),
                    FastMath.abs(expected * solver.getRelativeAccuracy()));
        result = solver.solve(100, f, min, max);
        Assert.assertEquals(expected, result, tolerance);

        min = -0.9; max = -0.2; expected = -0.5;
        tolerance = FastMath.max(solver.getAbsoluteAccuracy(),
                    FastMath.abs(expected * solver.getRelativeAccuracy()));
        result = solver.solve(100, f, min, max);
        Assert.assertEquals(expected, result, tolerance);
    }

    @Test
    public void testQuinticFunction() {
        UnivariateFunction f = new QuinticFunction();
        UnivariateSolver solver = new MullerSolver();
        double min, max, expected, result, tolerance;

        min = -0.4; max = 0.2; expected = 0.0;
        tolerance = FastMath.max(solver.getAbsoluteAccuracy(),
                    FastMath.abs(expected * solver.getRelativeAccuracy()));
        result = solver.solve(100, f, min, max);
        Assert.assertEquals(expected, result, tolerance);

        min = 0.75; max = 1.5; expected = 1.0;
        tolerance = FastMath.max(solver.getAbsoluteAccuracy(),
                    FastMath.abs(expected * solver.getRelativeAccuracy()));
        result = solver.solve(100, f, min, max);
        Assert.assertEquals(expected, result, tolerance);

        min = -0.9; max = -0.2; expected = -0.5;
        tolerance = FastMath.max(solver.getAbsoluteAccuracy(),
                    FastMath.abs(expected * solver.getRelativeAccuracy()));
        result = solver.solve(100, f, min, max);
        Assert.assertEquals(expected, result, tolerance);
    }

    @Test
    public void testQuinticFunction() {
        UnivariateFunction f = new QuinticFunction();
        UnivariateSolver solver = new RiddersSolver();
        double min, max, expected, result, tolerance;

        min = -0.4; max = 0.2; expected = 0.0;
        tolerance = FastMath.max(solver.getAbsoluteAccuracy(),
                    FastMath.abs(expected * solver.getRelativeAccuracy()));
        result = solver.solve(100, f, min, max);
        Assert.assertEquals(expected, result, tolerance);

        min = 0.75; max = 1.5; expected = 1.0;
        tolerance = FastMath.max(solver.getAbsoluteAccuracy(),
                    FastMath.abs(expected * solver.getRelativeAccuracy()));
        result = solver.solve(100, f, min, max);
        Assert.assertEquals(expected, result, tolerance);

        min = -0.9; max = -0.2; expected = -0.5;
        tolerance = FastMath.max(solver.getAbsoluteAccuracy(),
                    FastMath.abs(expected * solver.getRelativeAccuracy()));
        result = solver.solve(100, f, min, max);
        Assert.assertEquals(expected, result, tolerance);
    }

[2150, 'src/test/java', 'org.apache.commons.math3.distribution', 'HypergeometricDistributionTest', 'testDegenerateNoFailures', 110, 125]

[2150, 'src/test/java', 'org.apache.commons.math3.distribution', 'HypergeometricDistributionTest', 'testDegenerateNoSuccesses', 128, 143]

[2150, 'src/test/java', 'org.apache.commons.math3.distribution', 'HypergeometricDistributionTest', 'testDegenerateFullSample', 146, 161]

    @Test
    public void testDegenerateNoFailures() {
        HypergeometricDistribution dist = new HypergeometricDistribution(5,5,3);
        setDistribution(dist);
        setCumulativeTestPoints(new int[] {-1, 0, 1, 3, 10 });
        setCumulativeTestValues(new double[] {0d, 0d, 0d, 1d, 1d});
        setDensityTestPoints(new int[] {-1, 0, 1, 3, 10});
        setDensityTestValues(new double[] {0d, 0d, 0d, 1d, 0d});
        setInverseCumulativeTestPoints(new double[] {0.1d, 0.5d});
        setInverseCumulativeTestValues(new int[] {3, 3});
        verifyDensities();
        verifyCumulativeProbabilities();
        verifyInverseCumulativeProbabilities();
        Assert.assertEquals(dist.getSupportLowerBound(), 3);
        Assert.assertEquals(dist.getSupportUpperBound(), 3);
    }

    @Test
    public void testDegenerateNoSuccesses() {
        HypergeometricDistribution dist = new HypergeometricDistribution(5,0,3);
        setDistribution(dist);
        setCumulativeTestPoints(new int[] {-1, 0, 1, 3, 10 });
        setCumulativeTestValues(new double[] {0d, 1d, 1d, 1d, 1d});
        setDensityTestPoints(new int[] {-1, 0, 1, 3, 10});
        setDensityTestValues(new double[] {0d, 1d, 0d, 0d, 0d});
        setInverseCumulativeTestPoints(new double[] {0.1d, 0.5d});
        setInverseCumulativeTestValues(new int[] {0, 0});
        verifyDensities();
        verifyCumulativeProbabilities();
        verifyInverseCumulativeProbabilities();
        Assert.assertEquals(dist.getSupportLowerBound(), 0);
        Assert.assertEquals(dist.getSupportUpperBound(), 0);
    }

    @Test
    public void testDegenerateFullSample() {
        HypergeometricDistribution dist = new HypergeometricDistribution(5,3,5);
        setDistribution(dist);
        setCumulativeTestPoints(new int[] {-1, 0, 1, 3, 10 });
        setCumulativeTestValues(new double[] {0d, 0d, 0d, 1d, 1d});
        setDensityTestPoints(new int[] {-1, 0, 1, 3, 10});
        setDensityTestValues(new double[] {0d, 0d, 0d, 1d, 0d});
        setInverseCumulativeTestPoints(new double[] {0.1d, 0.5d});
        setInverseCumulativeTestValues(new int[] {3, 3});
        verifyDensities();
        verifyCumulativeProbabilities();
        verifyInverseCumulativeProbabilities();
        Assert.assertEquals(dist.getSupportLowerBound(), 3);
        Assert.assertEquals(dist.getSupportUpperBound(), 3);
    }

[2203, 'src/test/java', 'org.apache.commons.math3.optimization.direct', 'SimplexOptimizerMultiDirectionalTest', 'testMinimize1', 30, 43]

[2203, 'src/test/java', 'org.apache.commons.math3.optimization.direct', 'SimplexOptimizerNelderMeadTest', 'testMinimize1', 35, 48]

[2203, 'src/test/java', 'org.apache.commons.math3.optimization.direct', 'SimplexOptimizerNelderMeadTest', 'testMaximize1', 65, 78]

    @Test
    public void testMinimize1() {
        SimplexOptimizer optimizer = new SimplexOptimizer(1e-11, 1e-30);
        optimizer.setSimplex(new MultiDirectionalSimplex(new double[] { 0.2, 0.2 }));
        final FourExtrema fourExtrema = new FourExtrema();

        final PointValuePair optimum
            = optimizer.optimize(200, fourExtrema, GoalType.MINIMIZE, new double[] { -3, 0 });
        Assert.assertEquals(fourExtrema.xM, optimum.getPoint()[0], 4e-6);
        Assert.assertEquals(fourExtrema.yP, optimum.getPoint()[1], 3e-6);
        Assert.assertEquals(fourExtrema.valueXmYp, optimum.getValue(), 8e-13);
        Assert.assertTrue(optimizer.getEvaluations() > 120);
        Assert.assertTrue(optimizer.getEvaluations() < 150);
    }

    @Test
    public void testMinimize1() {
        SimplexOptimizer optimizer = new SimplexOptimizer(1e-10, 1e-30);
        optimizer.setSimplex(new NelderMeadSimplex(new double[] { 0.2, 0.2 }));
        final FourExtrema fourExtrema = new FourExtrema();

        final PointValuePair optimum
            = optimizer.optimize(100, fourExtrema, GoalType.MINIMIZE, new double[] { -3, 0 });
        Assert.assertEquals(fourExtrema.xM, optimum.getPoint()[0], 2e-7);
        Assert.assertEquals(fourExtrema.yP, optimum.getPoint()[1], 2e-5);
        Assert.assertEquals(fourExtrema.valueXmYp, optimum.getValue(), 6e-12);
        Assert.assertTrue(optimizer.getEvaluations() > 60);
        Assert.assertTrue(optimizer.getEvaluations() < 90);
    }

    @Test
    public void testMaximize1() {
        SimplexOptimizer optimizer = new SimplexOptimizer(1e-10, 1e-30);
        optimizer.setSimplex(new NelderMeadSimplex(new double[] { 0.2, 0.2 }));
        final FourExtrema fourExtrema = new FourExtrema();

        final PointValuePair optimum
            = optimizer.optimize(100, fourExtrema, GoalType.MAXIMIZE, new double[] { -3, 0 });
        Assert.assertEquals(fourExtrema.xM, optimum.getPoint()[0], 1e-5);
        Assert.assertEquals(fourExtrema.yM, optimum.getPoint()[1], 3e-6);
        Assert.assertEquals(fourExtrema.valueXmYm, optimum.getValue(), 3e-12);
        Assert.assertTrue(optimizer.getEvaluations() > 60);
        Assert.assertTrue(optimizer.getEvaluations() < 90);
    }

[2232, 'src/test/java', 'org.apache.commons.math3.ode.nonstiff', 'DormandPrince54IntegratorTest', 'testMinStep', 53, 74]

[2232, 'src/test/java', 'org.apache.commons.math3.ode.nonstiff', 'DormandPrince853IntegratorTest', 'testMinStep', 129, 150]

[2232, 'src/test/java', 'org.apache.commons.math3.ode.nonstiff', 'HighamHall54IntegratorTest', 'testMinStep', 76, 97]

  @Test(expected=NumberIsTooSmallException.class)
  public void testMinStep()
      throws DimensionMismatchException, NumberIsTooSmallException,
             MaxCountExceededException, NoBracketingException {

      TestProblem1 pb = new TestProblem1();
      double minStep = 0.1 * (pb.getFinalTime() - pb.getInitialTime());
      double maxStep = pb.getFinalTime() - pb.getInitialTime();
      double[] vecAbsoluteTolerance = { 1.0e-15, 1.0e-16 };
      double[] vecRelativeTolerance = { 1.0e-15, 1.0e-16 };

      FirstOrderIntegrator integ = new DormandPrince54Integrator(minStep, maxStep,
                                                                 vecAbsoluteTolerance,
                                                                 vecRelativeTolerance);
      TestProblemHandler handler = new TestProblemHandler(pb, integ);
      integ.addStepHandler(handler);
      integ.integrate(pb,
                      pb.getInitialTime(), pb.getInitialState(),
                      pb.getFinalTime(), new double[pb.getDimension()]);
      Assert.fail("an exception should have been thrown");

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

      FirstOrderIntegrator integ = new DormandPrince853Integrator(minStep, maxStep,
                                                                  vecAbsoluteTolerance,
                                                                  vecRelativeTolerance);
      TestProblemHandler handler = new TestProblemHandler(pb, integ);
      integ.addStepHandler(handler);
      integ.integrate(pb,
                      pb.getInitialTime(), pb.getInitialState(),
                      pb.getFinalTime(), new double[pb.getDimension()]);
      Assert.fail("an exception should have been thrown");

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

      FirstOrderIntegrator integ = new HighamHall54Integrator(minStep, maxStep,
                                                              vecAbsoluteTolerance,
                                                              vecRelativeTolerance);
      TestProblemHandler handler = new TestProblemHandler(pb, integ);
      integ.addStepHandler(handler);
      integ.integrate(pb,
                      pb.getInitialTime(), pb.getInitialState(),
                      pb.getFinalTime(), new double[pb.getDimension()]);
      Assert.fail("an exception should have been thrown");

  }

[2266, 'src/test/java', 'org.apache.commons.math3.linear', 'RealVectorFormatAbstractTest', 'testSimpleNoDecimals', 45, 51]

[2266, 'src/test/java', 'org.apache.commons.math3.linear', 'RealVectorFormatAbstractTest', 'testNonDefaultSetting', 113, 119]

[2266, 'src/test/java', 'org.apache.commons.math3.linear', 'RealVectorFormatAbstractTest', 'testParseSimpleNoDecimals', 166, 172]

    @Test
    public void testSimpleNoDecimals() {
        ArrayRealVector c = new ArrayRealVector(new double[] {1, 1, 1});
        String expected = "{1; 1; 1}";
        String actual = realVectorFormat.format(c);
        Assert.assertEquals(expected, actual);
    }

    @Test
    public void testNonDefaultSetting() {
        ArrayRealVector c = new ArrayRealVector(new double[] {1, 1, 1});
        String expected = "[1 : 1 : 1]";
        String actual = realVectorFormatSquare.format(c);
        Assert.assertEquals(expected, actual);
    }

    @Test
    public void testParseSimpleNoDecimals() {
        String source = "{1; 1; 1}";
        ArrayRealVector expected = new ArrayRealVector(new double[] {1, 1, 1});
        ArrayRealVector actual = realVectorFormat.parse(source);
        Assert.assertEquals(expected, actual);
    }

[2267, 'src/test/java', 'org.apache.commons.math3.linear', 'SchurTransformerTest', 'testPOrthogonal', 69, 74]

[2267, 'src/test/java', 'org.apache.commons.math3.linear', 'SchurTransformerTest', 'testPTOrthogonal', 76, 81]

[2267, 'src/test/java', 'org.apache.commons.math3.linear', 'SchurTransformerTest', 'testSchurForm', 83, 88]

    @Test
    public void testPOrthogonal() {
        checkOrthogonal(new SchurTransformer(MatrixUtils.createRealMatrix(testSquare5)).getP());
        checkOrthogonal(new SchurTransformer(MatrixUtils.createRealMatrix(testSquare3)).getP());
        checkOrthogonal(new SchurTransformer(MatrixUtils.createRealMatrix(testRandom)).getP());
    }

    @Test
    public void testPTOrthogonal() {
        checkOrthogonal(new SchurTransformer(MatrixUtils.createRealMatrix(testSquare5)).getPT());
        checkOrthogonal(new SchurTransformer(MatrixUtils.createRealMatrix(testSquare3)).getPT());
        checkOrthogonal(new SchurTransformer(MatrixUtils.createRealMatrix(testRandom)).getPT());
    }

    @Test
    public void testSchurForm() {
        checkSchurForm(new SchurTransformer(MatrixUtils.createRealMatrix(testSquare5)).getT());
        checkSchurForm(new SchurTransformer(MatrixUtils.createRealMatrix(testSquare3)).getT());
        checkSchurForm(new SchurTransformer(MatrixUtils.createRealMatrix(testRandom)).getT());
    }

[2309, 'src/test/java', 'org.apache.commons.math3.complex', 'ComplexTest', 'testCosInf', 748, 758]

[2309, 'src/test/java', 'org.apache.commons.math3.complex', 'ComplexTest', 'testCoshInf', 772, 782]

[2309, 'src/test/java', 'org.apache.commons.math3.complex', 'ComplexTest', 'testSinInf', 1003, 1013]

[2309, 'src/test/java', 'org.apache.commons.math3.complex', 'ComplexTest', 'testSinhInf', 1032, 1042]

    @Test
    public void testCosInf() {
        TestUtils.assertSame(infNegInf, oneInf.cos());
        TestUtils.assertSame(infInf, oneNegInf.cos());
        TestUtils.assertSame(Complex.NaN, infOne.cos());
        TestUtils.assertSame(Complex.NaN, negInfOne.cos());
        TestUtils.assertSame(Complex.NaN, infInf.cos());
        TestUtils.assertSame(Complex.NaN, infNegInf.cos());
        TestUtils.assertSame(Complex.NaN, negInfInf.cos());
        TestUtils.assertSame(Complex.NaN, negInfNegInf.cos());
    }

    @Test
    public void testCoshInf() {
        TestUtils.assertSame(Complex.NaN, oneInf.cosh());
        TestUtils.assertSame(Complex.NaN, oneNegInf.cosh());
        TestUtils.assertSame(infInf, infOne.cosh());
        TestUtils.assertSame(infNegInf, negInfOne.cosh());
        TestUtils.assertSame(Complex.NaN, infInf.cosh());
        TestUtils.assertSame(Complex.NaN, infNegInf.cosh());
        TestUtils.assertSame(Complex.NaN, negInfInf.cosh());
        TestUtils.assertSame(Complex.NaN, negInfNegInf.cosh());
    }

    @Test
    public void testSinInf() {
        TestUtils.assertSame(infInf, oneInf.sin());
        TestUtils.assertSame(infNegInf, oneNegInf.sin());
        TestUtils.assertSame(Complex.NaN, infOne.sin());
        TestUtils.assertSame(Complex.NaN, negInfOne.sin());
        TestUtils.assertSame(Complex.NaN, infInf.sin());
        TestUtils.assertSame(Complex.NaN, infNegInf.sin());
        TestUtils.assertSame(Complex.NaN, negInfInf.sin());
        TestUtils.assertSame(Complex.NaN, negInfNegInf.sin());
    }

    @Test
    public void testSinhInf() {
        TestUtils.assertSame(Complex.NaN, oneInf.sinh());
        TestUtils.assertSame(Complex.NaN, oneNegInf.sinh());
        TestUtils.assertSame(infInf, infOne.sinh());
        TestUtils.assertSame(negInfInf, negInfOne.sinh());
        TestUtils.assertSame(Complex.NaN, infInf.sinh());
        TestUtils.assertSame(Complex.NaN, infNegInf.sinh());
        TestUtils.assertSame(Complex.NaN, negInfInf.sinh());
        TestUtils.assertSame(Complex.NaN, negInfNegInf.sinh());
    }

[2390, 'src/test/java', 'org.apache.commons.math3.linear', 'RealMatrixFormatAbstractTest', 'testForgottenPrefix', 340, 346]

[2390, 'src/test/java', 'org.apache.commons.math3.linear', 'RealMatrixFormatAbstractTest', 'testForgottenSeparator', 348, 354]

[2390, 'src/test/java', 'org.apache.commons.math3.linear', 'RealMatrixFormatAbstractTest', 'testForgottenSuffix', 356, 362]

    @Test
    public void testForgottenPrefix() {
        ParsePosition pos = new ParsePosition(0);
        final String source = "1; 1; 1]";
        Assert.assertNull("Should not parse <"+source+">", realMatrixFormat.parse(source, pos));
        Assert.assertEquals(0, pos.getErrorIndex());
    }

    @Test
    public void testForgottenSeparator() {
        ParsePosition pos = new ParsePosition(0);
        final String source = "{{1, 1 1}}";
        Assert.assertNull("Should not parse <"+source+">", realMatrixFormat.parse(source, pos));
        Assert.assertEquals(7, pos.getErrorIndex());
    }

    @Test
    public void testForgottenSuffix() {
        ParsePosition pos = new ParsePosition(0);
        final String source = "{{1, 1, 1 ";
        Assert.assertNull("Should not parse <"+source+">", realMatrixFormat.parse(source, pos));
        Assert.assertEquals(9, pos.getErrorIndex());
    }

[2391, 'src/test/java', 'org.apache.commons.math3.linear', 'BiDiagonalTransformerTest', 'testUOrthogonal', 77, 82]

[2391, 'src/test/java', 'org.apache.commons.math3.linear', 'BiDiagonalTransformerTest', 'testVOrthogonal', 84, 89]

[2391, 'src/test/java', 'org.apache.commons.math3.linear', 'BiDiagonalTransformerTest', 'testBBiDiagonal', 97, 102]

[2391, 'src/test/java', 'org.apache.commons.math3.linear', 'SingularValueDecompositionTest', 'testUOrthogonal', 148, 153]

[2391, 'src/test/java', 'org.apache.commons.math3.linear', 'SingularValueDecompositionTest', 'testVOrthogonal', 156, 161]

    @Test
    public void testUOrthogonal() {
        checkOrthogonal(new BiDiagonalTransformer(MatrixUtils.createRealMatrix(testSquare)).getU());
        checkOrthogonal(new BiDiagonalTransformer(MatrixUtils.createRealMatrix(testNonSquare)).getU());
        checkOrthogonal(new BiDiagonalTransformer(MatrixUtils.createRealMatrix(testNonSquare).transpose()).getU());
    }

    @Test
    public void testVOrthogonal() {
        checkOrthogonal(new BiDiagonalTransformer(MatrixUtils.createRealMatrix(testSquare)).getV());
        checkOrthogonal(new BiDiagonalTransformer(MatrixUtils.createRealMatrix(testNonSquare)).getV());
        checkOrthogonal(new BiDiagonalTransformer(MatrixUtils.createRealMatrix(testNonSquare).transpose()).getV());
    }

    @Test
    public void testBBiDiagonal() {
        checkBiDiagonal(new BiDiagonalTransformer(MatrixUtils.createRealMatrix(testSquare)).getB());
        checkBiDiagonal(new BiDiagonalTransformer(MatrixUtils.createRealMatrix(testNonSquare)).getB());
        checkBiDiagonal(new BiDiagonalTransformer(MatrixUtils.createRealMatrix(testNonSquare).transpose()).getB());
    }

    @Test
    public void testUOrthogonal() {
        checkOrthogonal(new SingularValueDecomposition(MatrixUtils.createRealMatrix(testSquare)).getU());
        checkOrthogonal(new SingularValueDecomposition(MatrixUtils.createRealMatrix(testNonSquare)).getU());
        checkOrthogonal(new SingularValueDecomposition(MatrixUtils.createRealMatrix(testNonSquare).transpose()).getU());
    }

    @Test
    public void testVOrthogonal() {
        checkOrthogonal(new SingularValueDecomposition(MatrixUtils.createRealMatrix(testSquare)).getV());
        checkOrthogonal(new SingularValueDecomposition(MatrixUtils.createRealMatrix(testNonSquare)).getV());
        checkOrthogonal(new SingularValueDecomposition(MatrixUtils.createRealMatrix(testNonSquare).transpose()).getV());
    }

[2394, 'src/test/java', 'org.apache.commons.math3.complex', 'ComplexTest', 'testScalarAddNaN', 139, 145]

[2394, 'src/test/java', 'org.apache.commons.math3.complex', 'ComplexTest', 'testScalarDivideNaN', 271, 277]

[2394, 'src/test/java', 'org.apache.commons.math3.complex', 'ComplexTest', 'testScalarMultiplyNaN', 400, 406]

[2394, 'src/test/java', 'org.apache.commons.math3.complex', 'ComplexTest', 'testScalarSubtractNaN', 473, 479]

    @Test
    public void testScalarAddNaN() {
        Complex x = new Complex(3.0, 4.0);
        double yDouble = Double.NaN;
        Complex yComplex = new Complex(yDouble);
        Assert.assertEquals(x.add(yComplex), x.add(yDouble));
    }

    @Test
    public void testScalarDivideNaN() {
        Complex x = new Complex(3.0, 4.0);
        double yDouble = Double.NaN;
        Complex yComplex = new Complex(yDouble);
        Assert.assertEquals(x.divide(yComplex), x.divide(yDouble));
    }

    @Test
    public void testScalarMultiplyNaN() {
        Complex x = new Complex(3.0, 4.0);
        double yDouble = Double.NaN;
        Complex yComplex = new Complex(yDouble);
        Assert.assertEquals(x.multiply(yComplex), x.multiply(yDouble));
    }

    @Test
    public void testScalarSubtractNaN() {
        Complex x = new Complex(3.0, 4.0);
        double yDouble = Double.NaN;
        Complex yComplex = new Complex(yDouble);
        Assert.assertEquals(x.subtract(yComplex), x.subtract(yDouble));
    }

[2396, 'src/test/java', 'org.apache.commons.math3.geometry.euclidean.twod', 'Vector2DFormatAbstractTest', 'testParseSimpleWithDecimals', 176, 185]

[2396, 'src/test/java', 'org.apache.commons.math3.geometry.euclidean.twod', 'Vector2DFormatAbstractTest', 'testParseSimpleWithDecimalsTrunc', 187, 196]

[2396, 'src/test/java', 'org.apache.commons.math3.geometry.euclidean.twod', 'Vector2DFormatAbstractTest', 'testParseNegativeZ', 220, 229]

[2396, 'src/test/java', 'org.apache.commons.math3.geometry.euclidean.twod', 'Vector2DFormatAbstractTest', 'testParseNonDefaultSetting', 253, 262]

    @Test
    public void testParseSimpleWithDecimals() throws MathParseException {
        String source =
            "{1" + getDecimalCharacter() +
            "23; 1" + getDecimalCharacter() +
            "43}";
        Vector2D expected = new Vector2D(1.23, 1.43);
        Vector2D actual = vector2DFormat.parse(source);
        Assert.assertEquals(expected, actual);
    }

    @Test
    public void testParseSimpleWithDecimalsTrunc() throws MathParseException {
        String source =
            "{1" + getDecimalCharacter() +
            "2323; 1" + getDecimalCharacter() +
            "4343}";
        Vector2D expected = new Vector2D(1.2323, 1.4343);
        Vector2D actual = vector2DFormat.parse(source);
        Assert.assertEquals(expected, actual);
    }

    @Test
    public void testParseNegativeZ() throws MathParseException {
        String source =
            "{1" + getDecimalCharacter() +
            "2323; 1" + getDecimalCharacter() +
            "4343}";
        Vector2D expected = new Vector2D(1.2323, 1.4343);
        Vector2D actual = vector2DFormat.parse(source);
        Assert.assertEquals(expected, actual);
    }

    @Test
    public void testParseNonDefaultSetting() throws MathParseException {
        String source =
            "[1" + getDecimalCharacter() +
            "2323 : 1" + getDecimalCharacter() +
            "4343]";
        Vector2D expected = new Vector2D(1.2323, 1.4343);
        Vector2D actual = vector2DFormatSquare.parse(source);
        Assert.assertEquals(expected, actual);
    }

[2459, 'src/test/java', 'org.apache.commons.math3.linear', 'BlockFieldMatrixTest', 'testOperate', 375, 387]

[2459, 'src/test/java', 'org.apache.commons.math3.linear', 'BlockFieldMatrixTest', 'testPremultiplyVector', 445, 458]

[2459, 'src/test/java', 'org.apache.commons.math3.linear', 'FieldMatrixImplTest', 'testOperate', 268, 280]

[2459, 'src/test/java', 'org.apache.commons.math3.linear', 'FieldMatrixImplTest', 'testPremultiplyVector', 308, 321]

    @Test
    public void testOperate() {
        FieldMatrix<Fraction> m = new BlockFieldMatrix<Fraction>(id);
        TestUtils.assertEquals(testVector, m.operate(testVector));
        TestUtils.assertEquals(testVector, m.operate(new ArrayFieldVector<Fraction>(testVector)).toArray());
        m = new BlockFieldMatrix<Fraction>(bigSingular);
        try {
            m.operate(testVector);
            Assert.fail("Expecting illegalArgumentException");
        } catch (MathIllegalArgumentException ex) {
            // ignored
        }
    }

    @Test
    public void testPremultiplyVector() {
        FieldMatrix<Fraction> m = new BlockFieldMatrix<Fraction>(testData);
        TestUtils.assertEquals(m.preMultiply(testVector), preMultTest);
        TestUtils.assertEquals(m.preMultiply(new ArrayFieldVector<Fraction>(testVector).getData()),
                               preMultTest);
        m = new BlockFieldMatrix<Fraction>(bigSingular);
        try {
            m.preMultiply(testVector);
            Assert.fail("expecting MathIllegalArgumentException");
        } catch (MathIllegalArgumentException ex) {
            // ignored
        }
    }

    @Test
    public void testOperate() {
        FieldMatrix<Fraction> m = new Array2DRowFieldMatrix<Fraction>(id);
        TestUtils.assertEquals(testVector, m.operate(testVector));
        TestUtils.assertEquals(testVector, m.operate(new ArrayFieldVector<Fraction>(testVector)).toArray());
        m = new Array2DRowFieldMatrix<Fraction>(bigSingular);
        try {
            m.operate(testVector);
            Assert.fail("Expecting illegalArgumentException");
        } catch (MathIllegalArgumentException ex) {
            // ignored
        }
    }

    @Test
    public void testPremultiplyVector() {
        FieldMatrix<Fraction> m = new Array2DRowFieldMatrix<Fraction>(testData);
        TestUtils.assertEquals(m.preMultiply(testVector), preMultTest);
        TestUtils.assertEquals(m.preMultiply(new ArrayFieldVector<Fraction>(testVector).getData()),
                               preMultTest);
        m = new Array2DRowFieldMatrix<Fraction>(bigSingular);
        try {
            m.preMultiply(testVector);
            Assert.fail("expecting MathIllegalArgumentException");
        } catch (MathIllegalArgumentException ex) {
            // ignored
        }
    }

[2619, 'src/test/java', 'org.apache.commons.math3.distribution', 'LogNormalDistributionTest', 'testGetScale', 165, 169]

[2619, 'src/test/java', 'org.apache.commons.math3.distribution', 'LogNormalDistributionTest', 'testGetShape', 171, 175]

[2619, 'src/test/java', 'org.apache.commons.math3.distribution', 'NormalDistributionTest', 'testGetMean', 123, 127]

[2619, 'src/test/java', 'org.apache.commons.math3.distribution', 'NormalDistributionTest', 'testGetStandardDeviation', 129, 133]

[2619, 'src/test/java', 'org.apache.commons.math3.distribution', 'ParetoDistributionTest', 'testGetScale', 136, 140]

[2619, 'src/test/java', 'org.apache.commons.math3.distribution', 'ParetoDistributionTest', 'testGetShape', 142, 146]

    @Test
    public void testGetScale() {
        LogNormalDistribution distribution = (LogNormalDistribution)getDistribution();
        Assert.assertEquals(2.1, distribution.getScale(), 0);
    }

    @Test
    public void testGetShape() {
        LogNormalDistribution distribution = (LogNormalDistribution)getDistribution();
        Assert.assertEquals(1.4, distribution.getShape(), 0);
    }

    @Test
    public void testGetMean() {
        NormalDistribution distribution = (NormalDistribution) getDistribution();
        Assert.assertEquals(2.1, distribution.getMean(), 0);
    }

    @Test
    public void testGetStandardDeviation() {
        NormalDistribution distribution = (NormalDistribution) getDistribution();
        Assert.assertEquals(1.4, distribution.getStandardDeviation(), 0);
    }

    @Test
    public void testGetScale() {
        ParetoDistribution distribution = (ParetoDistribution)getDistribution();
        Assert.assertEquals(2.1, distribution.getScale(), 0);
    }

    @Test
    public void testGetShape() {
        ParetoDistribution distribution = (ParetoDistribution)getDistribution();
        Assert.assertEquals(1.4, distribution.getShape(), 0);
    }

[2683, 'src/test/java', 'org.apache.commons.math3.fraction', 'BigFractionFormatTest', 'testNumeratorFormat', 277, 292]

[2683, 'src/test/java', 'org.apache.commons.math3.fraction', 'BigFractionFormatTest', 'testDenominatorFormat', 294, 309]

[2683, 'src/test/java', 'org.apache.commons.math3.fraction', 'FractionFormatTest', 'testNumeratorFormat', 296, 311]

[2683, 'src/test/java', 'org.apache.commons.math3.fraction', 'FractionFormatTest', 'testDenominatorFormat', 313, 328]

    @Test
    public void testNumeratorFormat() {
        NumberFormat old = properFormat.getNumeratorFormat();
        NumberFormat nf = NumberFormat.getInstance();
        nf.setParseIntegerOnly(true);
        properFormat.setNumeratorFormat(nf);
        Assert.assertEquals(nf, properFormat.getNumeratorFormat());
        properFormat.setNumeratorFormat(old);

        old = improperFormat.getNumeratorFormat();
        nf = NumberFormat.getInstance();
        nf.setParseIntegerOnly(true);
        improperFormat.setNumeratorFormat(nf);
        Assert.assertEquals(nf, improperFormat.getNumeratorFormat());
        improperFormat.setNumeratorFormat(old);
    }

    @Test
    public void testDenominatorFormat() {
        NumberFormat old = properFormat.getDenominatorFormat();
        NumberFormat nf = NumberFormat.getInstance();
        nf.setParseIntegerOnly(true);
        properFormat.setDenominatorFormat(nf);
        Assert.assertEquals(nf, properFormat.getDenominatorFormat());
        properFormat.setDenominatorFormat(old);

        old = improperFormat.getDenominatorFormat();
        nf = NumberFormat.getInstance();
        nf.setParseIntegerOnly(true);
        improperFormat.setDenominatorFormat(nf);
        Assert.assertEquals(nf, improperFormat.getDenominatorFormat());
        improperFormat.setDenominatorFormat(old);
    }

    @Test
    public void testNumeratorFormat() {
        NumberFormat old = properFormat.getNumeratorFormat();
        NumberFormat nf = NumberFormat.getInstance();
        nf.setParseIntegerOnly(true);
        properFormat.setNumeratorFormat(nf);
        Assert.assertEquals(nf, properFormat.getNumeratorFormat());
        properFormat.setNumeratorFormat(old);

        old = improperFormat.getNumeratorFormat();
        nf = NumberFormat.getInstance();
        nf.setParseIntegerOnly(true);
        improperFormat.setNumeratorFormat(nf);
        Assert.assertEquals(nf, improperFormat.getNumeratorFormat());
        improperFormat.setNumeratorFormat(old);
    }

    @Test
    public void testDenominatorFormat() {
        NumberFormat old = properFormat.getDenominatorFormat();
        NumberFormat nf = NumberFormat.getInstance();
        nf.setParseIntegerOnly(true);
        properFormat.setDenominatorFormat(nf);
        Assert.assertEquals(nf, properFormat.getDenominatorFormat());
        properFormat.setDenominatorFormat(old);

        old = improperFormat.getDenominatorFormat();
        nf = NumberFormat.getInstance();
        nf.setParseIntegerOnly(true);
        improperFormat.setDenominatorFormat(nf);
        Assert.assertEquals(nf, improperFormat.getDenominatorFormat());
        improperFormat.setDenominatorFormat(old);
    }

[2775, 'src/test/java', 'org.apache.commons.math3.fraction', 'BigFractionTest', 'testIntValue', 250, 257]

[2775, 'src/test/java', 'org.apache.commons.math3.fraction', 'BigFractionTest', 'testLongValue', 259, 266]

[2775, 'src/test/java', 'org.apache.commons.math3.fraction', 'FractionTest', 'testIntValue', 206, 213]

[2775, 'src/test/java', 'org.apache.commons.math3.fraction', 'FractionTest', 'testLongValue', 215, 222]

    @Test
    public void testIntValue() {
        BigFraction first = new BigFraction(1, 2);
        BigFraction second = new BigFraction(3, 2);

        Assert.assertEquals(0, first.intValue());
        Assert.assertEquals(1, second.intValue());
    }

    @Test
    public void testLongValue() {
        BigFraction first = new BigFraction(1, 2);
        BigFraction second = new BigFraction(3, 2);

        Assert.assertEquals(0L, first.longValue());
        Assert.assertEquals(1L, second.longValue());
    }

    @Test
    public void testIntValue() {
        Fraction first = new Fraction(1, 2);
        Fraction second = new Fraction(3, 2);

        Assert.assertEquals(0, first.intValue());
        Assert.assertEquals(1, second.intValue());
    }

    @Test
    public void testLongValue() {
        Fraction first = new Fraction(1, 2);
        Fraction second = new Fraction(3, 2);

        Assert.assertEquals(0L, first.longValue());
        Assert.assertEquals(1L, second.longValue());
    }

[2779, 'src/test/java', 'org.apache.commons.math3.complex', 'ComplexFormatAbstractTest', 'testNegativeBoth', 109, 115]

[2779, 'src/test/java', 'org.apache.commons.math3.complex', 'ComplexFormatAbstractTest', 'testParseNegativeBoth', 218, 224]

[2779, 'src/test/java', 'org.apache.commons.math3.complex', 'ComplexFormatAbstractTest', 'testParseDifferentImaginaryChar', 242, 248]

    @Test
    public void testNegativeBoth() {
        Complex c = new Complex(-1.232323232323, -1.434343434343);
        String expected = "-1" + getDecimalCharacter() + "2323232323 - 1" + getDecimalCharacter() + "4343434343i";
        String actual = complexFormat.format(c);
        Assert.assertEquals(expected, actual);
    }

    @Test
    public void testParseNegativeBoth() {
        String source = "-1" + getDecimalCharacter() + "232323232323 - 1" + getDecimalCharacter() + "434343434343i";
        Complex expected = new Complex(-1.232323232323, -1.434343434343);
        Complex actual = complexFormat.parse(source);
        Assert.assertEquals(expected, actual);
    }

    @Test
    public void testParseDifferentImaginaryChar() {
        String source = "-1" + getDecimalCharacter() + "2323 - 1" + getDecimalCharacter() + "4343j";
        Complex expected = new Complex(-1.2323, -1.4343);
        Complex actual = complexFormatJ.parse(source);
        Assert.assertEquals(expected, actual);
    }

[2872, 'src/test/java', 'org.apache.commons.math3.analysis.differentiation', 'SparseGradientTest', 'testLogExp', 572, 580]

[2872, 'src/test/java', 'org.apache.commons.math3.analysis.differentiation', 'SparseGradientTest', 'testLog1pExpm1', 582, 590]

[2872, 'src/test/java', 'org.apache.commons.math3.analysis.differentiation', 'SparseGradientTest', 'testSinAsin', 615, 623]

[2872, 'src/test/java', 'org.apache.commons.math3.analysis.differentiation', 'SparseGradientTest', 'testCosAcos', 625, 633]

[2872, 'src/test/java', 'org.apache.commons.math3.analysis.differentiation', 'SparseGradientTest', 'testTanAtan', 635, 643]

[2872, 'src/test/java', 'org.apache.commons.math3.analysis.differentiation', 'SparseGradientTest', 'testSinhAsinh', 734, 742]

[2872, 'src/test/java', 'org.apache.commons.math3.analysis.differentiation', 'SparseGradientTest', 'testCoshAcosh', 744, 752]

[2872, 'src/test/java', 'org.apache.commons.math3.analysis.differentiation', 'SparseGradientTest', 'testTanhAtanh', 754, 762]

    @Test
    public void testLogExp() {
        for (double x = 0.1; x < 1.2; x += 0.001) {
            SparseGradient sgX = SparseGradient.createVariable(0, x);
            SparseGradient rebuiltX = sgX.exp().log();
            SparseGradient zero = rebuiltX.subtract(sgX);
            checkF0F1(zero, 0.0, 0.0);
        }
    }

    @Test
    public void testLog1pExpm1() {
        for (double x = 0.1; x < 1.2; x += 0.001) {
            SparseGradient sgX = SparseGradient.createVariable(0, x);
            SparseGradient rebuiltX = sgX.expm1().log1p();
            SparseGradient zero = rebuiltX.subtract(sgX);
            checkF0F1(zero, 0.0, 0.0);
        }
    }

    @Test
    public void testSinAsin() {
        for (double x = 0.1; x < 1.2; x += 0.001) {
            SparseGradient sgX = SparseGradient.createVariable(0, x);
            SparseGradient rebuiltX = sgX.sin().asin();
            SparseGradient zero = rebuiltX.subtract(sgX);
            checkF0F1(zero, 0.0, 0.0);
        }
    }

    @Test
    public void testCosAcos() {
        for (double x = 0.1; x < 1.2; x += 0.001) {
            SparseGradient sgX = SparseGradient.createVariable(0, x);
            SparseGradient rebuiltX = sgX.cos().acos();
            SparseGradient zero = rebuiltX.subtract(sgX);
            checkF0F1(zero, 0.0, 0.0);
        }
    }

    @Test
    public void testTanAtan() {
        for (double x = 0.1; x < 1.2; x += 0.001) {
            SparseGradient sgX = SparseGradient.createVariable(0, x);
            SparseGradient rebuiltX = sgX.tan().atan();
            SparseGradient zero = rebuiltX.subtract(sgX);
            checkF0F1(zero, 0.0, 0.0);
        }
    }

    @Test
    public void testSinhAsinh() {
        for (double x = 0.1; x < 1.2; x += 0.001) {
            SparseGradient sgX = SparseGradient.createVariable(0, x);
            SparseGradient rebuiltX = sgX.sinh().asinh();
            SparseGradient zero = rebuiltX.subtract(sgX);
            checkF0F1(zero, 0.0, 0.0);
        }
    }

    @Test
    public void testCoshAcosh() {
        for (double x = 0.1; x < 1.2; x += 0.001) {
            SparseGradient sgX = SparseGradient.createVariable(0, x);
            SparseGradient rebuiltX = sgX.cosh().acosh();
            SparseGradient zero = rebuiltX.subtract(sgX);
            checkF0F1(zero, 0.0, 0.0);
        }
    }

    @Test
    public void testTanhAtanh() {
        for (double x = 0.1; x < 1.2; x += 0.001) {
            SparseGradient sgX = SparseGradient.createVariable(0, x);
            SparseGradient rebuiltX = sgX.tanh().atanh();
            SparseGradient zero = rebuiltX.subtract(sgX);
            checkF0F1(zero, 0.0, 0.0);
        }
    }

[3002, 'src/test/java', 'org.apache.commons.math3.geometry.euclidean.threed', 'Vector3DFormatAbstractTest', 'testNegativeX', 78, 88]

[3002, 'src/test/java', 'org.apache.commons.math3.geometry.euclidean.threed', 'Vector3DFormatAbstractTest', 'testNegativeY', 90, 100]

[3002, 'src/test/java', 'org.apache.commons.math3.geometry.euclidean.threed', 'Vector3DFormatAbstractTest', 'testNegativeZ', 102, 112]

    @Test
    public void testNegativeX() {
        Vector3D c = new Vector3D(-1.232323232323, 1.43, 1.63);
        String expected =
            "{-1"    + getDecimalCharacter() +
            "2323232323; 1" + getDecimalCharacter() +
            "43; 1" + getDecimalCharacter() +
            "63}";
        String actual = vector3DFormat.format(c);
        Assert.assertEquals(expected, actual);
    }

    @Test
    public void testNegativeY() {
        Vector3D c = new Vector3D(1.23, -1.434343434343, 1.63);
        String expected =
            "{1"    + getDecimalCharacter() +
            "23; -1" + getDecimalCharacter() +
            "4343434343; 1" + getDecimalCharacter() +
            "63}";
        String actual = vector3DFormat.format(c);
        Assert.assertEquals(expected, actual);
    }

    @Test
    public void testNegativeZ() {
        Vector3D c = new Vector3D(1.23, 1.43, -1.633333333333);
        String expected =
            "{1"    + getDecimalCharacter() +
            "23; 1" + getDecimalCharacter() +
            "43; -1" + getDecimalCharacter() +
            "6333333333}";
        String actual = vector3DFormat.format(c);
        Assert.assertEquals(expected, actual);
    }

[3034, 'src/test/java', 'org.apache.commons.math3.analysis.interpolation', 'LoessInterpolatorTest', 'testNotAllFiniteReal1', 182, 185]

[3034, 'src/test/java', 'org.apache.commons.math3.analysis.interpolation', 'LoessInterpolatorTest', 'testNotAllFiniteReal2', 187, 190]

[3034, 'src/test/java', 'org.apache.commons.math3.analysis.interpolation', 'LoessInterpolatorTest', 'testNotAllFiniteReal3', 192, 195]

[3034, 'src/test/java', 'org.apache.commons.math3.analysis.interpolation', 'LoessInterpolatorTest', 'testNotAllFiniteReal4', 197, 200]

[3034, 'src/test/java', 'org.apache.commons.math3.analysis.interpolation', 'LoessInterpolatorTest', 'testNotAllFiniteReal5', 202, 205]

[3034, 'src/test/java', 'org.apache.commons.math3.analysis.interpolation', 'LoessInterpolatorTest', 'testNotAllFiniteReal6', 207, 210]

    @Test(expected=NotFiniteNumberException.class)
    public void testNotAllFiniteReal1() {
        new LoessInterpolator().smooth(new double[] {1,2,Double.NaN}, new double[] {3,4,5});
    }

    @Test(expected=NotFiniteNumberException.class)
    public void testNotAllFiniteReal2() {
        new LoessInterpolator().smooth(new double[] {1,2,Double.POSITIVE_INFINITY}, new double[] {3,4,5});
    }

    @Test(expected=NotFiniteNumberException.class)
    public void testNotAllFiniteReal3() {
        new LoessInterpolator().smooth(new double[] {1,2,Double.NEGATIVE_INFINITY}, new double[] {3,4,5});
    }

    @Test(expected=NotFiniteNumberException.class)
    public void testNotAllFiniteReal4() {
        new LoessInterpolator().smooth(new double[] {3,4,5}, new double[] {1,2,Double.NaN});
    }

    @Test(expected=NotFiniteNumberException.class)
    public void testNotAllFiniteReal5() {
        new LoessInterpolator().smooth(new double[] {3,4,5}, new double[] {1,2,Double.POSITIVE_INFINITY});
    }

    @Test(expected=NotFiniteNumberException.class)
    public void testNotAllFiniteReal6() {
        new LoessInterpolator().smooth(new double[] {3,4,5}, new double[] {1,2,Double.NEGATIVE_INFINITY});
    }

[3063, 'src/test/java', 'org.apache.commons.math3.random', 'RandomDataGeneratorTest', 'testNextBinomial', 1110, 1134]

[3063, 'src/test/java', 'org.apache.commons.math3.random', 'RandomDataGeneratorTest', 'testNextPascal', 1162, 1185]

[3063, 'src/test/java', 'org.apache.commons.math3.random', 'RandomDataGeneratorTest', 'testNextZipf', 1187, 1210]

    @Test
    public void testNextBinomial() {
        BinomialDistributionTest testInstance = new BinomialDistributionTest();
        int[] densityPoints = testInstance.makeDensityTestPoints();
        double[] densityValues = testInstance.makeDensityTestValues();
        int sampleSize = 1000;
        int length = TestUtils.eliminateZeroMassPoints(densityPoints, densityValues);
        BinomialDistribution distribution = (BinomialDistribution) testInstance.makeDistribution();
        double[] expectedCounts = new double[length];
        long[] observedCounts = new long[length];
        for (int i = 0; i < length; i++) {
            expectedCounts[i] = sampleSize * densityValues[i];
        }
        randomData.reSeed(1000);
        for (int i = 0; i < sampleSize; i++) {
          int value = randomData.nextBinomial(distribution.getNumberOfTrials(),
                  distribution.getProbabilityOfSuccess());
          for (int j = 0; j < length; j++) {
              if (value == densityPoints[j]) {
                  observedCounts[j]++;
              }
          }
        }
        TestUtils.assertChiSquareAccept(densityPoints, expectedCounts, observedCounts, .001);
    }

    @Test
    public void testNextPascal() {
        PascalDistributionTest testInstance = new PascalDistributionTest();
        int[] densityPoints = testInstance.makeDensityTestPoints();
        double[] densityValues = testInstance.makeDensityTestValues();
        int sampleSize = 1000;
        int length = TestUtils.eliminateZeroMassPoints(densityPoints, densityValues);
        PascalDistribution distribution = (PascalDistribution) testInstance.makeDistribution();
        double[] expectedCounts = new double[length];
        long[] observedCounts = new long[length];
        for (int i = 0; i < length; i++) {
            expectedCounts[i] = sampleSize * densityValues[i];
        }
        randomData.reSeed(1000);
        for (int i = 0; i < sampleSize; i++) {
          int value = randomData.nextPascal(distribution.getNumberOfSuccesses(), distribution.getProbabilityOfSuccess());
          for (int j = 0; j < length; j++) {
              if (value == densityPoints[j]) {
                  observedCounts[j]++;
              }
          }
        }
        TestUtils.assertChiSquareAccept(densityPoints, expectedCounts, observedCounts, .001);
    }

    @Test
    public void testNextZipf() {
        ZipfDistributionTest testInstance = new ZipfDistributionTest();
        int[] densityPoints = testInstance.makeDensityTestPoints();
        double[] densityValues = testInstance.makeDensityTestValues();
        int sampleSize = 1000;
        int length = TestUtils.eliminateZeroMassPoints(densityPoints, densityValues);
        ZipfDistribution distribution = (ZipfDistribution) testInstance.makeDistribution();
        double[] expectedCounts = new double[length];
        long[] observedCounts = new long[length];
        for (int i = 0; i < length; i++) {
            expectedCounts[i] = sampleSize * densityValues[i];
        }
        randomData.reSeed(1000);
        for (int i = 0; i < sampleSize; i++) {
          int value = randomData.nextZipf(distribution.getNumberOfElements(), distribution.getExponent());
          for (int j = 0; j < length; j++) {
              if (value == densityPoints[j]) {
                  observedCounts[j]++;
              }
          }
        }
        TestUtils.assertChiSquareAccept(densityPoints, expectedCounts, observedCounts, .001);
    }

[3067, 'src/test/java', 'org.apache.commons.math3.ode.nonstiff', 'DormandPrince54IntegratorTest', 'testEvents', 198, 231]

[3067, 'src/test/java', 'org.apache.commons.math3.ode.nonstiff', 'DormandPrince853IntegratorTest', 'testEvents', 243, 275]

[3067, 'src/test/java', 'org.apache.commons.math3.ode.nonstiff', 'GraggBulirschStoerIntegratorTest', 'testEvents', 194, 226]

[3067, 'src/test/java', 'org.apache.commons.math3.ode.nonstiff', 'HighamHall54IntegratorTest', 'testEvents', 161, 194]

  @Test
  public void testEvents()
      throws DimensionMismatchException, NumberIsTooSmallException,
             MaxCountExceededException, NoBracketingException {

    TestProblem4 pb = new TestProblem4();
    double minStep = 0;
    double maxStep = pb.getFinalTime() - pb.getInitialTime();
    double scalAbsoluteTolerance = 1.0e-8;
    double scalRelativeTolerance = 0.01 * scalAbsoluteTolerance;

    FirstOrderIntegrator integ = new DormandPrince54Integrator(minStep, maxStep,
                                                               scalAbsoluteTolerance,
                                                               scalRelativeTolerance);
    TestProblemHandler handler = new TestProblemHandler(pb, integ);
    integ.addStepHandler(handler);
    EventHandler[] functions = pb.getEventsHandlers();
    double convergence = 1.0e-8 * maxStep;
    for (int l = 0; l < functions.length; ++l) {
      integ.addEventHandler(functions[l],
                                 Double.POSITIVE_INFINITY, convergence, 1000);
    }
    Assert.assertEquals(functions.length, integ.getEventHandlers().size());
    integ.integrate(pb,
                    pb.getInitialTime(), pb.getInitialState(),
                    pb.getFinalTime(), new double[pb.getDimension()]);

    Assert.assertTrue(handler.getMaximalValueError() < 5.0e-6);
    Assert.assertEquals(0, handler.getMaximalTimeError(), convergence);
    Assert.assertEquals(12.0, handler.getLastTime(), convergence);
    integ.clearEventHandlers();
    Assert.assertEquals(0, integ.getEventHandlers().size());

  }

  @Test
  public void testEvents()
      throws DimensionMismatchException, NumberIsTooSmallException,
             MaxCountExceededException, NoBracketingException {

    TestProblem4 pb = new TestProblem4();
    double minStep = 0;
    double maxStep = pb.getFinalTime() - pb.getInitialTime();
    double scalAbsoluteTolerance = 1.0e-9;
    double scalRelativeTolerance = 0.01 * scalAbsoluteTolerance;

    FirstOrderIntegrator integ = new DormandPrince853Integrator(minStep, maxStep,
                                                                scalAbsoluteTolerance,
                                                                scalRelativeTolerance);
    TestProblemHandler handler = new TestProblemHandler(pb, integ);
    integ.addStepHandler(handler);
    EventHandler[] functions = pb.getEventsHandlers();
    double convergence = 1.0e-8 * maxStep;
    for (int l = 0; l < functions.length; ++l) {
      integ.addEventHandler(functions[l], Double.POSITIVE_INFINITY, convergence, 1000);
    }
    Assert.assertEquals(functions.length, integ.getEventHandlers().size());
    integ.integrate(pb,
                    pb.getInitialTime(), pb.getInitialState(),
                    pb.getFinalTime(), new double[pb.getDimension()]);

    Assert.assertEquals(0, handler.getMaximalValueError(), 2.1e-7);
    Assert.assertEquals(0, handler.getMaximalTimeError(), convergence);
    Assert.assertEquals(12.0, handler.getLastTime(), convergence);
    integ.clearEventHandlers();
    Assert.assertEquals(0, integ.getEventHandlers().size());

  }

  @Test
  public void testEvents()
      throws DimensionMismatchException, NumberIsTooSmallException,
             MaxCountExceededException, NoBracketingException {

    TestProblem4 pb = new TestProblem4();
    double minStep = 0;
    double maxStep = pb.getFinalTime() - pb.getInitialTime();
    double scalAbsoluteTolerance = 1.0e-10;
    double scalRelativeTolerance = 0.01 * scalAbsoluteTolerance;

    FirstOrderIntegrator integ = new GraggBulirschStoerIntegrator(minStep, maxStep,
                                                                  scalAbsoluteTolerance,
                                                                  scalRelativeTolerance);
    TestProblemHandler handler = new TestProblemHandler(pb, integ);
    integ.addStepHandler(handler);
    EventHandler[] functions = pb.getEventsHandlers();
    double convergence = 1.0e-8 * maxStep;
    for (int l = 0; l < functions.length; ++l) {
      integ.addEventHandler(functions[l], Double.POSITIVE_INFINITY, convergence, 1000);
    }
    Assert.assertEquals(functions.length, integ.getEventHandlers().size());
    integ.integrate(pb,
                    pb.getInitialTime(), pb.getInitialState(),
                    pb.getFinalTime(), new double[pb.getDimension()]);

    Assert.assertTrue(handler.getMaximalValueError() < 4.0e-7);
    Assert.assertEquals(0, handler.getMaximalTimeError(), convergence);
    Assert.assertEquals(12.0, handler.getLastTime(), convergence);
    integ.clearEventHandlers();
    Assert.assertEquals(0, integ.getEventHandlers().size());

  }

  @Test
  public void testEvents()
      throws DimensionMismatchException, NumberIsTooSmallException,
             MaxCountExceededException, NoBracketingException {

    TestProblem4 pb = new TestProblem4();
    double minStep = 0;
    double maxStep = pb.getFinalTime() - pb.getInitialTime();
    double scalAbsoluteTolerance = 1.0e-8;
    double scalRelativeTolerance = 0.01 * scalAbsoluteTolerance;

    FirstOrderIntegrator integ = new HighamHall54Integrator(minStep, maxStep,
                                                            scalAbsoluteTolerance,
                                                            scalRelativeTolerance);
    TestProblemHandler handler = new TestProblemHandler(pb, integ);
    integ.addStepHandler(handler);
    EventHandler[] functions = pb.getEventsHandlers();
    double convergence = 1.0e-8 * maxStep;
    for (int l = 0; l < functions.length; ++l) {
      integ.addEventHandler(functions[l],
                                 Double.POSITIVE_INFINITY, convergence, 1000);
    }
    Assert.assertEquals(functions.length, integ.getEventHandlers().size());
    integ.integrate(pb,
                    pb.getInitialTime(), pb.getInitialState(),
                    pb.getFinalTime(), new double[pb.getDimension()]);

    Assert.assertTrue(handler.getMaximalValueError() < 1.0e-7);
    Assert.assertEquals(0, handler.getMaximalTimeError(), convergence);
    Assert.assertEquals(12.0, handler.getLastTime(), convergence);
    integ.clearEventHandlers();
    Assert.assertEquals(0, integ.getEventHandlers().size());

  }

[3184, 'src/test/java', 'org.apache.commons.math3.linear', 'DiagonalMatrixTest', 'testSetNonDiagonalEntry', 275, 279]

[3184, 'src/test/java', 'org.apache.commons.math3.linear', 'DiagonalMatrixTest', 'testAddNonDiagonalEntry', 288, 292]

[3184, 'src/test/java', 'org.apache.commons.math3.linear', 'DiagonalMatrixTest', 'testSetEntryOutOfRange', 315, 319]

    @Test(expected=NumberIsTooLargeException.class)
    public void testSetNonDiagonalEntry() {
        final DiagonalMatrix diag = new DiagonalMatrix(3);
        diag.setEntry(1, 2, 3.4);
    }

    @Test(expected=NumberIsTooLargeException.class)
    public void testAddNonDiagonalEntry() {
        final DiagonalMatrix diag = new DiagonalMatrix(3);
        diag.addToEntry(1, 2, 3.4);
    }

    @Test(expected=OutOfRangeException.class)
    public void testSetEntryOutOfRange() {
        final DiagonalMatrix diag = new DiagonalMatrix(3);
        diag.setEntry(3, 3, 3.4);
    }

[3253, 'src/test/java', 'org.apache.commons.math3.ode.nonstiff', 'DormandPrince54StepInterpolatorTest', 'serialization', 58, 108]

[3253, 'src/test/java', 'org.apache.commons.math3.ode.nonstiff', 'DormandPrince853StepInterpolatorTest', 'serialization', 58, 108]

[3253, 'src/test/java', 'org.apache.commons.math3.ode.nonstiff', 'HighamHall54StepInterpolatorTest', 'serialization', 58, 108]

  @Test
  public void serialization()
    throws IOException, ClassNotFoundException,
           DimensionMismatchException, NumberIsTooSmallException,
           MaxCountExceededException, NoBracketingException  {

    TestProblem3 pb = new TestProblem3(0.9);
    double minStep = 0;
    double maxStep = pb.getFinalTime() - pb.getInitialTime();
    double scalAbsoluteTolerance = 1.0e-8;
    double scalRelativeTolerance = scalAbsoluteTolerance;
    DormandPrince54Integrator integ = new DormandPrince54Integrator(minStep, maxStep,
                                                                    scalAbsoluteTolerance,
                                                                    scalRelativeTolerance);
    integ.addStepHandler(new ContinuousOutputModel());
    integ.integrate(pb,
                    pb.getInitialTime(), pb.getInitialState(),
                    pb.getFinalTime(), new double[pb.getDimension()]);

    ByteArrayOutputStream bos = new ByteArrayOutputStream();
    ObjectOutputStream    oos = new ObjectOutputStream(bos);
    for (StepHandler handler : integ.getStepHandlers()) {
        oos.writeObject(handler);
    }

    Assert.assertTrue(bos.size () > 135000);
    Assert.assertTrue(bos.size () < 145000);

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

    Assert.assertTrue(maxError < 7.0e-10);

  }

  @Test
  public void serialization()
    throws IOException, ClassNotFoundException,
           DimensionMismatchException, NumberIsTooSmallException,
           MaxCountExceededException, NoBracketingException {

    TestProblem3 pb = new TestProblem3(0.9);
    double minStep = 0;
    double maxStep = pb.getFinalTime() - pb.getInitialTime();
    double scalAbsoluteTolerance = 1.0e-8;
    double scalRelativeTolerance = scalAbsoluteTolerance;
    DormandPrince853Integrator integ = new DormandPrince853Integrator(minStep, maxStep,
                                                                      scalAbsoluteTolerance,
                                                                      scalRelativeTolerance);
    integ.addStepHandler(new ContinuousOutputModel());
    integ.integrate(pb,
                    pb.getInitialTime(), pb.getInitialState(),
                    pb.getFinalTime(), new double[pb.getDimension()]);

    ByteArrayOutputStream bos = new ByteArrayOutputStream();
    ObjectOutputStream    oos = new ObjectOutputStream(bos);
    for (StepHandler handler : integ.getStepHandlers()) {
        oos.writeObject(handler);
    }

    Assert.assertTrue(bos.size () > 90000);
    Assert.assertTrue(bos.size () < 100000);

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

    Assert.assertTrue(maxError < 2.4e-10);

  }

  @Test
  public void serialization()
    throws IOException, ClassNotFoundException,
           DimensionMismatchException, NumberIsTooSmallException,
           MaxCountExceededException, NoBracketingException {

    TestProblem3 pb = new TestProblem3(0.9);
    double minStep = 0;
    double maxStep = pb.getFinalTime() - pb.getInitialTime();
    double scalAbsoluteTolerance = 1.0e-8;
    double scalRelativeTolerance = scalAbsoluteTolerance;
    HighamHall54Integrator integ = new HighamHall54Integrator(minStep, maxStep,
                                                              scalAbsoluteTolerance,
                                                              scalRelativeTolerance);
    integ.addStepHandler(new ContinuousOutputModel());
    integ.integrate(pb,
                    pb.getInitialTime(), pb.getInitialState(),
                    pb.getFinalTime(), new double[pb.getDimension()]);

    ByteArrayOutputStream bos = new ByteArrayOutputStream();
    ObjectOutputStream    oos = new ObjectOutputStream(bos);
    for (StepHandler handler : integ.getStepHandlers()) {
        oos.writeObject(handler);
    }

    Assert.assertTrue(bos.size () > 185000);
    Assert.assertTrue(bos.size () < 195000);

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

    Assert.assertTrue(maxError < 1.6e-10);

  }

[3336, 'src/test/java', 'org.apache.commons.math3.ode.nonstiff', 'AdamsBashforthIntegratorTest', 'dimensionCheck', 46, 54]

[3336, 'src/test/java', 'org.apache.commons.math3.ode.nonstiff', 'AdamsMoultonIntegratorTest', 'dimensionCheck', 44, 54]

[3336, 'src/test/java', 'org.apache.commons.math3.ode.nonstiff', 'DormandPrince54IntegratorTest', 'testDimensionCheck', 41, 51]

[3336, 'src/test/java', 'org.apache.commons.math3.ode.nonstiff', 'GraggBulirschStoerIntegratorTest', 'testDimensionCheck', 42, 52]

    @Test(expected=DimensionMismatchException.class)
    public void dimensionCheck() throws NumberIsTooSmallException, DimensionMismatchException, MaxCountExceededException, NoBracketingException {
        TestProblem1 pb = new TestProblem1();
        FirstOrderIntegrator integ =
            new AdamsBashforthIntegrator(2, 0.0, 1.0, 1.0e-10, 1.0e-10);
        integ.integrate(pb,
                        0.0, new double[pb.getDimension()+10],
                        1.0, new double[pb.getDimension()+10]);
    }

    @Test(expected=DimensionMismatchException.class)
    public void dimensionCheck()
        throws DimensionMismatchException, NumberIsTooSmallException,
               MaxCountExceededException, NoBracketingException {
        TestProblem1 pb = new TestProblem1();
        FirstOrderIntegrator integ =
            new AdamsMoultonIntegrator(2, 0.0, 1.0, 1.0e-10, 1.0e-10);
        integ.integrate(pb,
                        0.0, new double[pb.getDimension()+10],
                        1.0, new double[pb.getDimension()+10]);
    }

  @Test(expected=DimensionMismatchException.class)
  public void testDimensionCheck()
      throws DimensionMismatchException, NumberIsTooSmallException,
             MaxCountExceededException, NoBracketingException {
      TestProblem1 pb = new TestProblem1();
      DormandPrince54Integrator integrator = new DormandPrince54Integrator(0.0, 1.0,
                                                                           1.0e-10, 1.0e-10);
      integrator.integrate(pb,
                           0.0, new double[pb.getDimension()+10],
                           1.0, new double[pb.getDimension()+10]);
  }

  @Test(expected=DimensionMismatchException.class)
  public void testDimensionCheck()
      throws DimensionMismatchException, NumberIsTooSmallException,
             MaxCountExceededException, NoBracketingException {
      TestProblem1 pb = new TestProblem1();
      AdaptiveStepsizeIntegrator integrator =
        new GraggBulirschStoerIntegrator(0.0, 1.0, 1.0e-10, 1.0e-10);
      integrator.integrate(pb,
                           0.0, new double[pb.getDimension()+10],
                           1.0, new double[pb.getDimension()+10]);
  }

[3340, 'src/test/java', 'org.apache.commons.math3.geometry.euclidean.oned', 'Vector1DTest', 'testDistance', 141, 148]

[3340, 'src/test/java', 'org.apache.commons.math3.geometry.euclidean.oned', 'Vector1DTest', 'testDistanceSq', 150, 158]

[3340, 'src/test/java', 'org.apache.commons.math3.geometry.euclidean.oned', 'Vector1DTest', 'testDistanceInf', 160, 167]

    @Test
    public void testDistance() {
        Vector1D v1 = new Vector1D(1);
        Vector1D v2 = new Vector1D(-4);
        Assert.assertEquals(0.0, Vector1D.distance(new Vector1D(-1), new Vector1D(-1)), 0);
        Assert.assertEquals(5.0, Vector1D.distance(v1, v2), 1.0e-12);
        Assert.assertEquals(v1.subtract(v2).getNorm(), Vector1D.distance(v1, v2), 1.0e-12);
    }

    @Test
    public void testDistanceSq() {
        Vector1D v1 = new Vector1D(1);
        Vector1D v2 = new Vector1D(-4);
        Assert.assertEquals(0.0, Vector1D.distanceSq(new Vector1D(-1), new Vector1D(-1)), 0);
        Assert.assertEquals(25.0, Vector1D.distanceSq(v1, v2), 1.0e-12);
        Assert.assertEquals(Vector1D.distance(v1, v2) * Vector1D.distance(v1, v2),
                            Vector1D.distanceSq(v1, v2), 1.0e-12);
  }

    @Test
    public void testDistanceInf() {
        Vector1D v1 = new Vector1D(1);
        Vector1D v2 = new Vector1D(-4);
        Assert.assertEquals(0.0, Vector1D.distanceInf(new Vector1D(-1), new Vector1D(-1)), 0);
        Assert.assertEquals(5.0, Vector1D.distanceInf(v1, v2), 1.0e-12);
        Assert.assertEquals(v1.subtract(v2).getNormInf(), Vector1D.distanceInf(v1, v2), 1.0e-12);
    }

[3415, 'src/test/java', 'org.apache.commons.math3.util', 'MathArraysTest', 'testL1DistanceDouble', 133, 138]

[3415, 'src/test/java', 'org.apache.commons.math3.util', 'MathArraysTest', 'testL2DistanceDouble', 147, 152]

[3415, 'src/test/java', 'org.apache.commons.math3.util', 'MathArraysTest', 'testLInfDistanceDouble', 161, 166]

    @Test
    public void testL1DistanceDouble() {
        double[] p1 = { 2.5,  0.0 };
        double[] p2 = { -0.5, 4.0 };
        Assert.assertTrue(Precision.equals(7.0, MathArrays.distance1(p1, p2), 1));
    }

    @Test
    public void testL2DistanceDouble() {
        double[] p1 = { 2.5,  0.0 };
        double[] p2 = { -0.5, 4.0 };
        Assert.assertTrue(Precision.equals(5.0, MathArrays.distance(p1, p2), 1));
    }

    @Test
    public void testLInfDistanceDouble() {
        double[] p1 = { 2.5,  0.0 };
        double[] p2 = { -0.5, 4.0 };
        Assert.assertTrue(Precision.equals(4.0, MathArrays.distanceInf(p1, p2), 1));
    }

[3434, 'src/test/java', 'org.apache.commons.math3.linear', 'SparseFieldMatrixTest', 'testGetRowMatrix', 477, 496]

[3434, 'src/test/java', 'org.apache.commons.math3.linear', 'SparseFieldMatrixTest', 'testGetColumnMatrix', 498, 517]

[3434, 'src/test/java', 'org.apache.commons.math3.linear', 'SparseFieldMatrixTest', 'testGetColumnVector', 540, 559]

    @Test
    public void testGetRowMatrix() {
        FieldMatrix<Fraction> m = createSparseMatrix(subTestData);
        FieldMatrix<Fraction> mRow0 = createSparseMatrix(subRow0);
        FieldMatrix<Fraction> mRow3 = createSparseMatrix(subRow3);
        Assert.assertEquals("Row0", mRow0, m.getRowMatrix(0));
        Assert.assertEquals("Row3", mRow3, m.getRowMatrix(3));
        try {
            m.getRowMatrix(-1);
            Assert.fail("Expecting OutOfRangeException");
        } catch (OutOfRangeException ex) {
            // expected
        }
        try {
            m.getRowMatrix(4);
            Assert.fail("Expecting OutOfRangeException");
        } catch (OutOfRangeException ex) {
            // expected
        }
    }

    @Test
    public void testGetColumnMatrix() {
        FieldMatrix<Fraction> m = createSparseMatrix(subTestData);
        FieldMatrix<Fraction> mColumn1 = createSparseMatrix(subColumn1);
        FieldMatrix<Fraction> mColumn3 = createSparseMatrix(subColumn3);
        Assert.assertEquals("Column1", mColumn1, m.getColumnMatrix(1));
        Assert.assertEquals("Column3", mColumn3, m.getColumnMatrix(3));
        try {
            m.getColumnMatrix(-1);
            Assert.fail("Expecting OutOfRangeException");
        } catch (OutOfRangeException ex) {
            // expected
        }
        try {
            m.getColumnMatrix(4);
            Assert.fail("Expecting OutOfRangeException");
        } catch (OutOfRangeException ex) {
            // expected
        }
    }

    @Test
    public void testGetColumnVector() {
        FieldMatrix<Fraction> m = createSparseMatrix(subTestData);
        FieldVector<Fraction> mColumn1 = columnToVector(subColumn1);
        FieldVector<Fraction> mColumn3 = columnToVector(subColumn3);
        Assert.assertEquals("Column1", mColumn1, m.getColumnVector(1));
        Assert.assertEquals("Column3", mColumn3, m.getColumnVector(3));
        try {
            m.getColumnVector(-1);
            Assert.fail("Expecting OutOfRangeException");
        } catch (OutOfRangeException ex) {
            // expected
        }
        try {
            m.getColumnVector(4);
            Assert.fail("Expecting OutOfRangeException");
        } catch (OutOfRangeException ex) {
            // expected
        }
    }

[3473, 'src/test/java', 'org.apache.commons.math3.linear', 'RealMatrixFormatAbstractTest', 'testParseSimpleWithDecimals', 214, 224]

[3473, 'src/test/java', 'org.apache.commons.math3.linear', 'RealMatrixFormatAbstractTest', 'testParseSimpleWithDecimalsTrunc', 226, 236]

[3473, 'src/test/java', 'org.apache.commons.math3.linear', 'RealMatrixFormatAbstractTest', 'testParseNonDefaultSetting', 274, 284]

    @Test
    public void testParseSimpleWithDecimals() {
        String source =
            "{{1" + getDecimalCharacter() +
            "23,1" + getDecimalCharacter() +
            "43,1" + getDecimalCharacter() +
            "63}}";
        RealMatrix expected = MatrixUtils.createRealMatrix(new double[][] {{1.23, 1.43, 1.63}});
        RealMatrix actual = realMatrixFormat.parse(source);
        Assert.assertEquals(expected, actual);
    }

    @Test
    public void testParseSimpleWithDecimalsTrunc() {
        String source =
            "{{1" + getDecimalCharacter() +
            "2323,1" + getDecimalCharacter() +
            "4343,1" + getDecimalCharacter() +
            "6333}}";
        RealMatrix expected = MatrixUtils.createRealMatrix(new double[][] {{1.2323, 1.4343, 1.6333}});
        RealMatrix actual = realMatrixFormat.parse(source);
        Assert.assertEquals(expected, actual);
    }

    @Test
    public void testParseNonDefaultSetting() {
        String source =
            "[1" + getDecimalCharacter() +
            "2323, 1" + getDecimalCharacter() +
            "4343, 1" + getDecimalCharacter() +
            "6333]";
        RealMatrix expected = MatrixUtils.createRealMatrix(new double[][] {{1.2323, 1.4343, 1.6333}});
        RealMatrix actual = realMatrixFormatOctave.parse(source);
        Assert.assertEquals(expected, actual);
    }

[3484, 'src/test/java', 'org.apache.commons.math3.complex', 'RootsOfUnityTest', 'testMathIllegalState1', 32, 36]

[3484, 'src/test/java', 'org.apache.commons.math3.complex', 'RootsOfUnityTest', 'testMathIllegalState2', 38, 42]

[3484, 'src/test/java', 'org.apache.commons.math3.complex', 'RootsOfUnityTest', 'testZeroNumberOfRoots', 50, 54]

    @Test(expected = MathIllegalStateException.class)
    public void testMathIllegalState1() {
        final RootsOfUnity roots = new RootsOfUnity();
        roots.getReal(0);
    }

    @Test(expected = MathIllegalStateException.class)
    public void testMathIllegalState2() {
        final RootsOfUnity roots = new RootsOfUnity();
        roots.getImaginary(0);
    }

    @Test(expected = ZeroException.class)
    public void testZeroNumberOfRoots() {
        final RootsOfUnity roots = new RootsOfUnity();
        roots.computeRoots(0);
    }

[3527, 'src/test/java', 'org.apache.commons.math3.complex', 'ComplexTest', 'testAdd', 99, 106]

[3527, 'src/test/java', 'org.apache.commons.math3.complex', 'ComplexTest', 'testDivide', 181, 188]

[3527, 'src/test/java', 'org.apache.commons.math3.complex', 'ComplexTest', 'testMultiply', 341, 348]

[3527, 'src/test/java', 'org.apache.commons.math3.complex', 'ComplexTest', 'testSubtract', 434, 441]

    @Test
    public void testAdd() {
        Complex x = new Complex(3.0, 4.0);
        Complex y = new Complex(5.0, 6.0);
        Complex z = x.add(y);
        Assert.assertEquals(8.0, z.getReal(), 1.0e-5);
        Assert.assertEquals(10.0, z.getImaginary(), 1.0e-5);
    }

    @Test
    public void testDivide() {
        Complex x = new Complex(3.0, 4.0);
        Complex y = new Complex(5.0, 6.0);
        Complex z = x.divide(y);
        Assert.assertEquals(39.0 / 61.0, z.getReal(), 1.0e-5);
        Assert.assertEquals(2.0 / 61.0, z.getImaginary(), 1.0e-5);
    }

    @Test
    public void testMultiply() {
        Complex x = new Complex(3.0, 4.0);
        Complex y = new Complex(5.0, 6.0);
        Complex z = x.multiply(y);
        Assert.assertEquals(-9.0, z.getReal(), 1.0e-5);
        Assert.assertEquals(38.0, z.getImaginary(), 1.0e-5);
    }

    @Test
    public void testSubtract() {
        Complex x = new Complex(3.0, 4.0);
        Complex y = new Complex(5.0, 6.0);
        Complex z = x.subtract(y);
        Assert.assertEquals(-2.0, z.getReal(), 1.0e-5);
        Assert.assertEquals(-2.0, z.getImaginary(), 1.0e-5);
    }

[3548, 'src/test/java', 'org.apache.commons.math3.linear', 'SparseRealMatrixTest', 'testGetRowMatrix', 468, 487]

[3548, 'src/test/java', 'org.apache.commons.math3.linear', 'SparseRealMatrixTest', 'testGetColumnMatrix', 489, 508]

[3548, 'src/test/java', 'org.apache.commons.math3.linear', 'SparseRealMatrixTest', 'testGetColumnVector', 531, 550]

    @Test
    public void testGetRowMatrix() {
        RealMatrix m = createSparseMatrix(subTestData);
        RealMatrix mRow0 = createSparseMatrix(subRow0);
        RealMatrix mRow3 = createSparseMatrix(subRow3);
        Assert.assertEquals("Row0", mRow0, m.getRowMatrix(0));
        Assert.assertEquals("Row3", mRow3, m.getRowMatrix(3));
        try {
            m.getRowMatrix(-1);
            Assert.fail("Expecting OutOfRangeException");
        } catch (OutOfRangeException ex) {
            // expected
        }
        try {
            m.getRowMatrix(4);
            Assert.fail("Expecting OutOfRangeException");
        } catch (OutOfRangeException ex) {
            // expected
        }
    }

    @Test
    public void testGetColumnMatrix() {
        RealMatrix m = createSparseMatrix(subTestData);
        RealMatrix mColumn1 = createSparseMatrix(subColumn1);
        RealMatrix mColumn3 = createSparseMatrix(subColumn3);
        Assert.assertEquals("Column1", mColumn1, m.getColumnMatrix(1));
        Assert.assertEquals("Column3", mColumn3, m.getColumnMatrix(3));
        try {
            m.getColumnMatrix(-1);
            Assert.fail("Expecting OutOfRangeException");
        } catch (OutOfRangeException ex) {
            // expected
        }
        try {
            m.getColumnMatrix(4);
            Assert.fail("Expecting OutOfRangeException");
        } catch (OutOfRangeException ex) {
            // expected
        }
    }

    @Test
    public void testGetColumnVector() {
        RealMatrix m = createSparseMatrix(subTestData);
        RealVector mColumn1 = columnToVector(subColumn1);
        RealVector mColumn3 = columnToVector(subColumn3);
        Assert.assertEquals("Column1", mColumn1, m.getColumnVector(1));
        Assert.assertEquals("Column3", mColumn3, m.getColumnVector(3));
        try {
            m.getColumnVector(-1);
            Assert.fail("Expecting OutOfRangeException");
        } catch (OutOfRangeException ex) {
            // expected
        }
        try {
            m.getColumnVector(4);
            Assert.fail("Expecting OutOfRangeException");
        } catch (OutOfRangeException ex) {
            // expected
        }
    }

[3549, 'src/test/java', 'org.apache.commons.math3.linear', 'CholeskySolverTest', 'testSolveDimensionErrors', 36, 59]

[3549, 'src/test/java', 'org.apache.commons.math3.linear', 'LUSolverTest', 'testSolveDimensionErrors', 74, 97]

[3549, 'src/test/java', 'org.apache.commons.math3.linear', 'LUSolverTest', 'testSolveSingularityErrors', 100, 123]

[3549, 'src/test/java', 'org.apache.commons.math3.linear', 'SingularValueSolverTest', 'testSolveDimensionErrors', 41, 64]

    @Test
    public void testSolveDimensionErrors() {
        DecompositionSolver solver =
            new CholeskyDecomposition(MatrixUtils.createRealMatrix(testData)).getSolver();
        RealMatrix b = MatrixUtils.createRealMatrix(new double[2][2]);
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
        try {
            solver.solve(new ArrayRealVectorTest.RealVectorTestImpl(b.getColumn(0)));
            Assert.fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException iae) {
            // expected behavior
        }
    }

    @Test
    public void testSolveDimensionErrors() {
        DecompositionSolver solver =
            new LUDecomposition(MatrixUtils.createRealMatrix(testData)).getSolver();
        RealMatrix b = MatrixUtils.createRealMatrix(new double[2][2]);
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
        try {
            solver.solve(new ArrayRealVectorTest.RealVectorTestImpl(b.getColumn(0)));
            Assert.fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException iae) {
            // expected behavior
        }
    }

    @Test
    public void testSolveSingularityErrors() {
        DecompositionSolver solver =
            new LUDecomposition(MatrixUtils.createRealMatrix(singular)).getSolver();
        RealMatrix b = MatrixUtils.createRealMatrix(new double[2][2]);
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
        try {
            solver.solve(new ArrayRealVectorTest.RealVectorTestImpl(b.getColumn(0)));
            Assert.fail("an exception should have been thrown");
        } catch (SingularMatrixException ime) {
            // expected behavior
        }
    }

    @Test
    public void testSolveDimensionErrors() {
        DecompositionSolver solver =
            new SingularValueDecomposition(MatrixUtils.createRealMatrix(testSquare)).getSolver();
        RealMatrix b = MatrixUtils.createRealMatrix(new double[3][2]);
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
        try {
            solver.solve(new ArrayRealVectorTest.RealVectorTestImpl(b.getColumn(0)));
            Assert.fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException iae) {
            // expected behavior
        }
    }

[3550, 'src/test/java', 'org.apache.commons.math3.linear', 'Array2DRowRealMatrixTest', 'testGetRowMatrix', 614, 635]

[3550, 'src/test/java', 'org.apache.commons.math3.linear', 'Array2DRowRealMatrixTest', 'testGetColumnMatrix', 658, 679]

[3550, 'src/test/java', 'org.apache.commons.math3.linear', 'BlockRealMatrixTest', 'testGetRowMatrix', 660, 679]

    @Test
    public void testGetRowMatrix() {
        RealMatrix m = new Array2DRowRealMatrix(subTestData);
        RealMatrix mRow0 = new Array2DRowRealMatrix(subRow0);
        RealMatrix mRow3 = new Array2DRowRealMatrix(subRow3);
        Assert.assertEquals("Row0", mRow0,
                m.getRowMatrix(0));
        Assert.assertEquals("Row3", mRow3,
                m.getRowMatrix(3));
        try {
            m.getRowMatrix(-1);
            Assert.fail("Expecting OutOfRangeException");
        } catch (OutOfRangeException ex) {
            // expected
        }
        try {
            m.getRowMatrix(4);
            Assert.fail("Expecting OutOfRangeException");
        } catch (OutOfRangeException ex) {
            // expected
        }
    }

    @Test
    public void testGetColumnMatrix() {
        RealMatrix m = new Array2DRowRealMatrix(subTestData);
        RealMatrix mColumn1 = new Array2DRowRealMatrix(subColumn1);
        RealMatrix mColumn3 = new Array2DRowRealMatrix(subColumn3);
        Assert.assertEquals("Column1", mColumn1,
                m.getColumnMatrix(1));
        Assert.assertEquals("Column3", mColumn3,
                m.getColumnMatrix(3));
        try {
            m.getColumnMatrix(-1);
            Assert.fail("Expecting OutOfRangeException");
        } catch (OutOfRangeException ex) {
            // expected
        }
        try {
            m.getColumnMatrix(4);
            Assert.fail("Expecting OutOfRangeException");
        } catch (OutOfRangeException ex) {
            // expected
        }
    }

    @Test
    public void testGetRowMatrix() {
        RealMatrix m     = new BlockRealMatrix(subTestData);
        RealMatrix mRow0 = new BlockRealMatrix(subRow0);
        RealMatrix mRow3 = new BlockRealMatrix(subRow3);
        Assert.assertEquals("Row0", mRow0, m.getRowMatrix(0));
        Assert.assertEquals("Row3", mRow3, m.getRowMatrix(3));
        try {
            m.getRowMatrix(-1);
            Assert.fail("Expecting OutOfRangeException");
        } catch (OutOfRangeException ex) {
            // expected
        }
        try {
            m.getRowMatrix(4);
            Assert.fail("Expecting OutOfRangeException");
        } catch (OutOfRangeException ex) {
            // expected
        }
    }

[3561, 'src/test/java', 'org.apache.commons.math3.geometry.euclidean.threed', 'Vector3DFormatAbstractTest', 'testParseSimpleWithDecimals', 184, 194]

[3561, 'src/test/java', 'org.apache.commons.math3.geometry.euclidean.threed', 'Vector3DFormatAbstractTest', 'testParseSimpleWithDecimalsTrunc', 196, 206]

[3561, 'src/test/java', 'org.apache.commons.math3.geometry.euclidean.threed', 'Vector3DFormatAbstractTest', 'testParseNonDefaultSetting', 268, 278]

    @Test
    public void testParseSimpleWithDecimals() throws MathParseException {
        String source =
            "{1" + getDecimalCharacter() +
            "23; 1" + getDecimalCharacter() +
            "43; 1" + getDecimalCharacter() +
            "63}";
        Vector3D expected = new Vector3D(1.23, 1.43, 1.63);
        Vector3D actual = vector3DFormat.parse(source);
        Assert.assertEquals(expected, actual);
    }

    @Test
    public void testParseSimpleWithDecimalsTrunc() throws MathParseException {
        String source =
            "{1" + getDecimalCharacter() +
            "2323; 1" + getDecimalCharacter() +
            "4343; 1" + getDecimalCharacter() +
            "6333}";
        Vector3D expected = new Vector3D(1.2323, 1.4343, 1.6333);
        Vector3D actual = vector3DFormat.parse(source);
        Assert.assertEquals(expected, actual);
    }

    @Test
    public void testParseNonDefaultSetting() throws MathParseException {
        String source =
            "[1" + getDecimalCharacter() +
            "2323 : 1" + getDecimalCharacter() +
            "4343 : 1" + getDecimalCharacter() +
            "6333]";
        Vector3D expected = new Vector3D(1.2323, 1.4343, 1.6333);
        Vector3D actual = vector3DFormatSquare.parse(source);
        Assert.assertEquals(expected, actual);
    }

[3590, 'src/test/java', 'org.apache.commons.math3.geometry.euclidean.twod', 'SubLineTest', 'testIntersectionInsideOutside', 115, 121]

[3590, 'src/test/java', 'org.apache.commons.math3.geometry.euclidean.twod', 'SubLineTest', 'testIntersectionBoundaryOutside', 131, 137]

[3590, 'src/test/java', 'org.apache.commons.math3.geometry.euclidean.twod', 'SubLineTest', 'testIntersectionOutsideOutside', 139, 145]

[3590, 'src/test/java', 'org.apache.commons.math3.geometry.euclidean.twod', 'SubLineTest', 'testIntersectionParallel', 147, 153]

    @Test
    public void testIntersectionInsideOutside() {
        SubLine sub1 = new SubLine(new Vector2D(1, 1), new Vector2D(3, 1), 1.0e-10);
        SubLine sub2 = new SubLine(new Vector2D(2, 0), new Vector2D(2, 0.5), 1.0e-10);
        Assert.assertNull(sub1.intersection(sub2, true));
        Assert.assertNull(sub1.intersection(sub2, false));
    }

    @Test
    public void testIntersectionBoundaryOutside() {
        SubLine sub1 = new SubLine(new Vector2D(1, 1), new Vector2D(2, 1), 1.0e-10);
        SubLine sub2 = new SubLine(new Vector2D(2, 0), new Vector2D(2, 0.5), 1.0e-10);
        Assert.assertNull(sub1.intersection(sub2, true));
        Assert.assertNull(sub1.intersection(sub2, false));
    }

    @Test
    public void testIntersectionOutsideOutside() {
        SubLine sub1 = new SubLine(new Vector2D(1, 1), new Vector2D(1.5, 1), 1.0e-10);
        SubLine sub2 = new SubLine(new Vector2D(2, 0), new Vector2D(2, 0.5), 1.0e-10);
        Assert.assertNull(sub1.intersection(sub2, true));
        Assert.assertNull(sub1.intersection(sub2, false));
    }

    @Test
    public void testIntersectionParallel() {
        final SubLine sub1 = new SubLine(new Vector2D(0, 1), new Vector2D(0, 2), 1.0e-10);
        final SubLine sub2 = new SubLine(new Vector2D(66, 3), new Vector2D(66, 4), 1.0e-10);
        Assert.assertNull(sub1.intersection(sub2, true));
        Assert.assertNull(sub1.intersection(sub2, false));
    }

[3645, 'src/test/java', 'org.apache.commons.math3.linear', 'RealMatrixFormatAbstractTest', 'testSimpleNoDecimals', 46, 52]

[3645, 'src/test/java', 'org.apache.commons.math3.linear', 'RealMatrixFormatAbstractTest', 'testNonDefaultSetting', 133, 139]

[3645, 'src/test/java', 'org.apache.commons.math3.linear', 'RealMatrixFormatAbstractTest', 'testParseSimpleNoDecimals', 184, 190]

    @Test
    public void testSimpleNoDecimals() {
        RealMatrix m = MatrixUtils.createRealMatrix(new double[][] {{1, 1, 1}, {1, 1, 1}});
        String expected = "{{1,1,1},{1,1,1}}";
        String actual = realMatrixFormat.format(m);
        Assert.assertEquals(expected, actual);
    }

    @Test
    public void testNonDefaultSetting() {
        RealMatrix m = MatrixUtils.createRealMatrix(new double[][] {{1, 1, 1}, {1, 1, 1}});
        String expected = "[1, 1, 1; 1, 1, 1]";
        String actual = realMatrixFormatOctave.format(m);
        Assert.assertEquals(expected, actual);
    }

    @Test
    public void testParseSimpleNoDecimals() {
        String source = "{{1, 1, 1}, {1, 1, 1}}";
        RealMatrix expected = MatrixUtils.createRealMatrix(new double[][] {{1, 1, 1}, {1, 1, 1}});
        RealMatrix actual = realMatrixFormat.parse(source);
        Assert.assertEquals(expected, actual);
    }

[3685, 'src/test/java', 'org.apache.commons.math3.ode.nonstiff', 'EulerIntegratorTest', 'testDimensionCheck', 44, 53]

[3685, 'src/test/java', 'org.apache.commons.math3.ode.nonstiff', 'GillIntegratorTest', 'testDimensionCheck', 44, 53]

[3685, 'src/test/java', 'org.apache.commons.math3.ode.nonstiff', 'MidpointIntegratorTest', 'testDimensionCheck', 44, 53]

[3685, 'src/test/java', 'org.apache.commons.math3.ode.nonstiff', 'ThreeEighthesIntegratorTest', 'testDimensionCheck', 44, 53]

  @Test(expected=DimensionMismatchException.class)
  public void testDimensionCheck()
      throws DimensionMismatchException, NumberIsTooSmallException,
             MaxCountExceededException, NoBracketingException {
      TestProblem1 pb = new TestProblem1();
      new EulerIntegrator(0.01).integrate(pb,
                                          0.0, new double[pb.getDimension()+10],
                                          1.0, new double[pb.getDimension()+10]);
        Assert.fail("an exception should have been thrown");
  }

  @Test(expected=DimensionMismatchException.class)
  public void testDimensionCheck()
      throws DimensionMismatchException, NumberIsTooSmallException,
             MaxCountExceededException, NoBracketingException {
      TestProblem1 pb = new TestProblem1();
      new GillIntegrator(0.01).integrate(pb,
                                         0.0, new double[pb.getDimension()+10],
                                         1.0, new double[pb.getDimension()+10]);
        Assert.fail("an exception should have been thrown");
  }

  @Test(expected=DimensionMismatchException.class)
  public void testDimensionCheck()
      throws DimensionMismatchException, NumberIsTooSmallException,
             MaxCountExceededException, NoBracketingException {
      TestProblem1 pb = new TestProblem1();
      new MidpointIntegrator(0.01).integrate(pb,
                                             0.0, new double[pb.getDimension()+10],
                                             1.0, new double[pb.getDimension()+10]);
        Assert.fail("an exception should have been thrown");
  }

  @Test(expected=DimensionMismatchException.class)
  public void testDimensionCheck()
      throws DimensionMismatchException, NumberIsTooSmallException,
             MaxCountExceededException, NoBracketingException {
      TestProblem1 pb = new TestProblem1();
      new ThreeEighthesIntegrator(0.01).integrate(pb,
                                                  0.0, new double[pb.getDimension()+10],
                                                  1.0, new double[pb.getDimension()+10]);
        Assert.fail("an exception should have been thrown");
  }

[3688, 'src/test/java', 'org.apache.commons.math3.linear', 'RealVectorFormatAbstractTest', 'testNegativeX', 77, 87]

[3688, 'src/test/java', 'org.apache.commons.math3.linear', 'RealVectorFormatAbstractTest', 'testNegativeY', 89, 99]

[3688, 'src/test/java', 'org.apache.commons.math3.linear', 'RealVectorFormatAbstractTest', 'testNegativeZ', 101, 111]

[3688, 'src/test/java', 'org.apache.commons.math3.linear', 'RealVectorFormatAbstractTest', 'testParseNegativeX', 211, 221]

[3688, 'src/test/java', 'org.apache.commons.math3.linear', 'RealVectorFormatAbstractTest', 'testParseNegativeY', 223, 233]

[3688, 'src/test/java', 'org.apache.commons.math3.linear', 'RealVectorFormatAbstractTest', 'testParseNegativeZ', 235, 245]

[3688, 'src/test/java', 'org.apache.commons.math3.linear', 'RealVectorFormatAbstractTest', 'testParseZeroX', 259, 269]

    @Test
    public void testNegativeX() {
        ArrayRealVector c = new ArrayRealVector(new double[] {-1.232323232323, 1.43, 1.63});
        String expected =
            "{-1"    + getDecimalCharacter() +
            "2323232323; 1" + getDecimalCharacter() +
            "43; 1" + getDecimalCharacter() +
            "63}";
        String actual = realVectorFormat.format(c);
        Assert.assertEquals(expected, actual);
    }

    @Test
    public void testNegativeY() {
        ArrayRealVector c = new ArrayRealVector(new double[] {1.23, -1.434343434343, 1.63});
        String expected =
            "{1"    + getDecimalCharacter() +
            "23; -1" + getDecimalCharacter() +
            "4343434343; 1" + getDecimalCharacter() +
            "63}";
        String actual = realVectorFormat.format(c);
        Assert.assertEquals(expected, actual);
    }

    @Test
    public void testNegativeZ() {
        ArrayRealVector c = new ArrayRealVector(new double[] {1.23, 1.43, -1.633333333333});
        String expected =
            "{1"    + getDecimalCharacter() +
            "23; 1" + getDecimalCharacter() +
            "43; -1" + getDecimalCharacter() +
            "6333333333}";
        String actual = realVectorFormat.format(c);
        Assert.assertEquals(expected, actual);
    }

    @Test
    public void testParseNegativeX() {
        String source =
            "{-1" + getDecimalCharacter() +
            "2323; 1" + getDecimalCharacter() +
            "4343; 1" + getDecimalCharacter() +
            "6333}";
        ArrayRealVector expected = new ArrayRealVector(new double[] {-1.2323, 1.4343, 1.6333});
        ArrayRealVector actual = realVectorFormat.parse(source);
        Assert.assertEquals(expected, actual);
    }

    @Test
    public void testParseNegativeY() {
        String source =
            "{1" + getDecimalCharacter() +
            "2323; -1" + getDecimalCharacter() +
            "4343; 1" + getDecimalCharacter() +
            "6333}";
        ArrayRealVector expected = new ArrayRealVector(new double[] {1.2323, -1.4343, 1.6333});
        ArrayRealVector actual = realVectorFormat.parse(source);
        Assert.assertEquals(expected, actual);
    }

    @Test
    public void testParseNegativeZ() {
        String source =
            "{1" + getDecimalCharacter() +
            "2323; 1" + getDecimalCharacter() +
            "4343; -1" + getDecimalCharacter() +
            "6333}";
        ArrayRealVector expected = new ArrayRealVector(new double[] {1.2323, 1.4343, -1.6333});
        ArrayRealVector actual = realVectorFormat.parse(source);
        Assert.assertEquals(expected, actual);
    }

    @Test
    public void testParseZeroX() {
        String source =
            "{0" + getDecimalCharacter() +
            "0; -1" + getDecimalCharacter() +
            "4343; 1" + getDecimalCharacter() +
            "6333}";
        ArrayRealVector expected = new ArrayRealVector(new double[] {0.0, -1.4343, 1.6333});
        ArrayRealVector actual = realVectorFormat.parse(source);
        Assert.assertEquals(expected, actual);
    }

