package ca.uwaterloo.jrefactoring.template;

import ca.uwaterloo.jrefactoring.node.RFNodeDifference;
import ca.uwaterloo.jrefactoring.utility.ASTNodeUtil;
import ca.uwaterloo.jrefactoring.utility.FileLogger;
import ca.uwaterloo.jrefactoring.utility.RenameUtil;
import ca.uwaterloo.jrefactoring.visitor.MethodVisitor;
import org.eclipse.core.runtime.IPath;
import org.eclipse.jdt.core.*;
import org.eclipse.jdt.core.dom.*;
import org.eclipse.jdt.core.dom.rewrite.ASTRewrite;
import org.eclipse.jdt.core.dom.rewrite.ListRewrite;
import org.eclipse.jface.text.Document;
import org.eclipse.text.edits.TextEdit;
import org.slf4j.Logger;

import java.io.File;
import java.util.*;

import static org.eclipse.jdt.core.dom.CompilationUnit.IMPORTS_PROPERTY;

public class RFTemplate {

    private static Logger log = FileLogger.getLogger(RFTemplate.class);
    private static final String TYPE_NAME_PREFIX = "T";
    private static final String CLAZZ_NAME_PREFIX = "clazz";
    private static final String CLASS_NAME = "Class";
    private static final String EXCEPTION_NAME = "Exception";
    private static final String THROWABLE_NAME = "Throwable";
    private static final String DEFAULT_LAMBDA_VARIABLE_NAME = "f";
    private static final String DEFAULT_LAMBDA_PARAMETER_NAME = "x";
    private static final String DEFAULT_LAMBDA_FUNCTION_CLASS_TEMPLATE = "java.util.function.To%sFunction";
    private static final String DEFAULT_ADAPTER_VARIABLE_NAME = "adapter";
    private static final String DEFAULT_ADAPTER_METHOD_NAME = "action";
    private static final String DEFAULT_TEMPLATE_CLASS_NAME = "TestTemplates";
    public static final String NULL_PARAMETER_OBJECT = "nullParameterObject";
    private static final String JAVA_LANG_CLASS = "java.lang.Class";
    private static final String ASSERT_KEYWORD = "assert";
    private static final String NEW_KEYWORD = "new";
    private static final String CASE_KEYWORD = "case";
    private static final String SHORT_KEYWORD = "short";
    private static final String FOR_KEYWORD = "for";
    private static final Set<String> keywords =
            new HashSet<>(Arrays.asList(ASSERT_KEYWORD, NEW_KEYWORD, CASE_KEYWORD, SHORT_KEYWORD, FOR_KEYWORD));
    private static final String KEYWORD_SUFFIX = "Action";
    private static final String OBJECT_NAME = "Object";
    private static final String JAVA_OBJECT_FULL_NAME = "java.lang.Object";
    private static final Set<String> javaObjectMethodNames = new HashSet<>(Arrays.asList("hashCode", "equals", "clone",
            "toString", "finalize"));
    private static final boolean ENABLE_LAMBDA = true;

    private AST ast;
    private MethodDeclaration templateMethod;
    private TypeDeclaration templateClass;
    private CompilationUnit templateCU;
    private TypeDeclaration adapter;
    private TypeDeclaration adapterImpl1;
    private TypeDeclaration adapterImpl2;
    private Map<TypePair, String> typeMap;
    private Map<String, TypePair> genericTypeMap;
    private Map<MethodInvocationPair, String> methodInvocationMap;
    private Map<String, String> clazzInstanceMap;
    private Map<String, String> nameMap1;
    private Map<String, String> nameMap2;
    private Map<String, Integer> parameterMap;
    private SingleVariableDeclaration lambdaVariable;
    private SingleVariableDeclaration adapterVariable;
    private Set<String> adapterTypes;
    private Map<ClassInstanceCreation, Type> instanceCreationTypeMap;
    private List<Expression> templateArguments1;
    private List<Expression> templateArguments2;
    private MethodDeclaration method1;
    private MethodDeclaration method2;
    private CompilationUnit compilationUnit1;
    private CompilationUnit compilationUnit2;
    private PackageDeclaration packageDeclaration1;
    private PackageDeclaration packageDeclaration2;
    private CompilationUnit adapterInterfaceCU;
    private CompilationUnit adapterImplCU1;
    private CompilationUnit adapterImplCU2;
    private List<Name> cuImports1;
    private List<Name> cuImports2;
    private List<Name> templateCUImports;
    private int clazzCount;
    private int typeCount;
    private int actionCount;
    private int variableCount;
    private boolean hasLambda;
    private boolean lambdaParameterAdded;
    private Type lambdaParameterType;
    private boolean hasAdapterVariable;
    private List<NodePair> unrefactoredList;
    private Map<String, Integer> adapterActionNameMap;
    private Map<String, Integer> genericTypeNameMap;
    private ICompilationUnit iCU1;
    private ICompilationUnit iCU2;
    private IPackageFragmentRoot packageFragmentRoot;
    private boolean refactorable;
    private boolean hasNullParameterObject;
    private Set<String> templateVariableNames;
    private boolean throwsAllException;
    private boolean accessIssue;

    public RFTemplate(AST ast, MethodDeclaration method1, MethodDeclaration method2,
                      String templateName, String adapterName, String[] adapterImplNamePair,
                      ICompilationUnit iCU1, ICompilationUnit iCU2) {

        assert adapterImplNamePair.length == 2;
        this.ast = ast;
        this.typeMap = new HashMap<>();
        this.genericTypeMap = new HashMap<>();
        this.methodInvocationMap = new HashMap<>();
        this.clazzInstanceMap = new HashMap<>();
        this.nameMap1 = new HashMap<>();
        this.nameMap2 = new HashMap<>();
        this.parameterMap = new HashMap<>();
        this.adapterTypes = new HashSet<>();
        this.instanceCreationTypeMap = new HashMap<>();
        this.clazzCount = 1;
        this.typeCount = 1;
        this.actionCount = 1;
        this.variableCount = 1;
        this.hasLambda = false;
        this.lambdaParameterAdded = false;
        this.lambdaParameterType = null;
        this.hasAdapterVariable = false;
        this.templateArguments1 = new ArrayList<>();
        this.templateArguments2 = new ArrayList<>();
        this.method1 = (MethodDeclaration) ASTNode.copySubtree(this.ast, method1);
        this.method2 = (MethodDeclaration) ASTNode.copySubtree(this.ast, method2);
        this.compilationUnit1 = (CompilationUnit) method1.getRoot();
        this.compilationUnit2 = (CompilationUnit) method2.getRoot();
        this.packageDeclaration1 = compilationUnit1.getPackage();
        this.packageDeclaration2 = compilationUnit2.getPackage();
        this.cuImports1 = new ArrayList<>();
        this.cuImports2 = new ArrayList<>();
        this.templateCUImports = new ArrayList<>();
        this.unrefactoredList = new ArrayList<>();
        this.adapterActionNameMap = new HashMap<>();
        this.genericTypeNameMap = new HashMap<>();
        this.iCU1 = iCU1;
        this.iCU2 = iCU2;
        this.packageFragmentRoot = (IPackageFragmentRoot) iCU1.getAncestor(3);
        this.refactorable = true;
        this.hasNullParameterObject = false;
        this.templateVariableNames = new HashSet<>();
        this.throwsAllException = false;
        this.accessIssue = false;
        init(templateName, adapterName, adapterImplNamePair);
    }

    private void init(String templateName, String adapterName, String[] adapterImplNamePair) {
        initTemplate(templateName);
        initAdapter(adapterName);
        initAdapterImpl(adapterImplNamePair[0], adapterImplNamePair[1]);
        initPackageDeclaration();
    }

    private void initTemplate(String templateName) {
        templateMethod = ast.newMethodDeclaration();
        templateMethod.setBody(ast.newBlock());
        templateMethod.setName(ast.newSimpleName(templateName));

        String name1 = compilationUnit1.getJavaElement().getElementName();
        String name2 = compilationUnit2.getJavaElement().getElementName();
        if (!packageDeclaration1.getName().getFullyQualifiedName().equals(packageDeclaration2.getName().getFullyQualifiedName())
                || !name1.equals(name2)) {
            templateMethod.modifiers().add(ast.newModifier(Modifier.ModifierKeyword.PUBLIC_KEYWORD));
            templateMethod.modifiers().add(ast.newModifier(Modifier.ModifierKeyword.STATIC_KEYWORD));

            templateClass = ast.newTypeDeclaration();
            templateClass.setName(ast.newSimpleName(templateName.substring(0, 1).toUpperCase() + templateName.substring(1)));
            templateClass.modifiers().add(ast.newModifier(Modifier.ModifierKeyword.PUBLIC_KEYWORD));
            templateClass.bodyDeclarations().add(templateMethod);

            templateCU = ast.newCompilationUnit();
            templateCU.types().add(templateClass);

        } else {
            templateClass = null;
            templateCU = null;
            templateMethod.modifiers().add(ast.newModifier(Modifier.ModifierKeyword.PUBLIC_KEYWORD));
        }

        if (method1.thrownExceptionTypes().size() > 0 || method2.thrownExceptionTypes().size() > 0) {
            throwsAllException = true;
        }

    }

    private void initAdapter(String adapterName) {
        // init lambda variable
        this.lambdaVariable = ast.newSingleVariableDeclaration();
        lambdaVariable.setName(ast.newSimpleName(DEFAULT_LAMBDA_VARIABLE_NAME));

        // init Adapter interface
        adapter = ast.newTypeDeclaration();
        Modifier publicModifier = ast.newModifier(Modifier.ModifierKeyword.PUBLIC_KEYWORD);
        adapter.modifiers().add(publicModifier);
        adapter.setInterface(true);
        adapter.setName(ast.newSimpleName(adapterName));

        // init Adapter CU
        adapterInterfaceCU = ast.newCompilationUnit();
        adapterInterfaceCU.types().add(adapter);

        // init adapter variable
        this.adapterVariable = ast.newSingleVariableDeclaration();
        Type type = ast.newSimpleType((SimpleName) ASTNode.copySubtree(ast, adapter.getName()));
        adapterVariable.setType(type);
        adapterVariable.setName(ast.newSimpleName(DEFAULT_ADAPTER_VARIABLE_NAME));
        //String adapterVariableName = adapterName.substring(0, 1).toLowerCase() + adapterName.substring(1);
        //adapterVariable.setName(ast.newSimpleName(adapterVariableName));

    }

    private void initAdapterImpl(String adapterImplName1, String adapterImplName2) {
        // init Adapter impl class
        adapterImpl1 = ast.newTypeDeclaration();
        adapterImpl2 = ast.newTypeDeclaration();
        Modifier publicModifier1 = ast.newModifier(Modifier.ModifierKeyword.PUBLIC_KEYWORD);
        Modifier publicModifier2 = ast.newModifier(Modifier.ModifierKeyword.PUBLIC_KEYWORD);
        adapterImpl1.modifiers().add(publicModifier1);
        adapterImpl2.modifiers().add(publicModifier2);
        adapterImpl1.setName(ast.newSimpleName(adapterImplName1));
        adapterImpl2.setName(ast.newSimpleName(adapterImplName2));

        Type interfaceType = ast.newSimpleType(ast.newSimpleName(adapter.getName().getIdentifier()));
        adapterImpl1.superInterfaceTypes().add(ASTNode.copySubtree(ast, interfaceType));
        adapterImpl2.superInterfaceTypes().add(ASTNode.copySubtree(ast, interfaceType));

        // init Adapter impl CU
        adapterImplCU1 = ast.newCompilationUnit();
        adapterImplCU2 = ast.newCompilationUnit();
        adapterImplCU1.types().add(adapterImpl1);
        adapterImplCU2.types().add(adapterImpl2);
    }

    private void initPackageDeclaration() {
        PackageDeclaration packageDeclaration = ast.newPackageDeclaration();
        String commonPackageName = ASTNodeUtil.getCommonPackageName(packageDeclaration1.getName().getFullyQualifiedName(),
                packageDeclaration2.getName().getFullyQualifiedName());
        if (commonPackageName != null && commonPackageName.length() > 0) {
            Name packageName = ASTNodeUtil.createPackageName(ast, commonPackageName);
            if (packageName != null) {
                packageDeclaration.setName(packageName);
                adapterInterfaceCU.setPackage(packageDeclaration);
                adapterImplCU1.setPackage((PackageDeclaration) ASTNode.copySubtree(ast, packageDeclaration));
                adapterImplCU2.setPackage((PackageDeclaration) ASTNode.copySubtree(ast, packageDeclaration));

                if (templateCU != null) {
                    templateCU.setPackage((PackageDeclaration) ASTNode.copySubtree(ast, packageDeclaration));
                    initImportDeclaration(commonPackageName);
                }
            }
        }
    }

    private void initImportDeclaration(String commonPackageName) {
        String templateName = commonPackageName + "." + templateClass.getName().getIdentifier();
        if (!commonPackageName.equals(packageDeclaration1.getName().getFullyQualifiedName())) {
            String adapterImplName1 = commonPackageName + "." + adapterImpl1.getName().getIdentifier();
            cuImports1.add(ASTNodeUtil.createPackageName(ast, templateName));
            cuImports1.add(ASTNodeUtil.createPackageName(ast, adapterImplName1));
        }
        if (!commonPackageName.equals(packageDeclaration2.getName().getFullyQualifiedName())) {
            String adapterImplName2 = commonPackageName + "." + adapterImpl2.getName().getIdentifier();
            cuImports2.add(ASTNodeUtil.createPackageName(ast, templateName));
            cuImports2.add(ASTNodeUtil.createPackageName(ast, adapterImplName2));
        }
    }

    public enum Pair {
        member1,
        member2
    }

    public AST getAst() {
        return ast;
    }

    public void setModifier(Modifier.ModifierKeyword modifier) {
        templateMethod.modifiers().add(ast.newModifier(modifier));
    }

    public boolean containsTypePair(TypePair typePair) {
        return typeMap.containsKey(typePair);
    }

    public boolean containsGenericType(String type) {
        return clazzInstanceMap.containsKey(type);
    }

    public String getClazz(String type) {
        return clazzInstanceMap.get(type);
    }

    public String getGenericTypeName(TypePair typePair) {
        return typeMap.get(typePair);
    }

    public ICompilationUnit getiCU1() {
        return iCU1;
    }

    public ICompilationUnit getiCU2() {
        return iCU2;
    }

    public boolean containsGenericNameInMap(String genericName) {
        return genericTypeMap.containsKey(genericName);
    }

    public TypePair getTypePairFromGenericMap(String genericName) {
        return genericTypeMap.get(genericName);
    }

    public boolean containTypeBound(String fromType, String toType) {
        List<TypeParameter> typeParameters = templateMethod.typeParameters();
        for (TypeParameter typeParameter : typeParameters) {
            if (typeParameter.getName().getIdentifier().equals(fromType)) {
                List<Type> typeBounds = typeParameter.typeBounds();
                for (Type typeBound : typeBounds) {
                    if (typeBound.isParameterizedType()) {
                        if (typeBound.toString().equals(toType)) {
                            return true;
                        }

                    } else {
                        if (((SimpleType) typeBound).getName().getFullyQualifiedName().equals(toType)) {
                            return true;
                        }
                    }
                }
            }
        }
        return false;
    }

    public List<Expression> getTemplateArguments1() {
        return templateArguments1;
    }

    public List<Expression> getTemplateArguments2() {
        return templateArguments2;
    }

    public CompilationUnit getTemplateCU() {
        return templateCU;
    }

    public boolean isRefactorable() {
        return refactorable && !hasUnrefactorableNodePair();
    }

    public void markAsUnrefactorable() {
        refactorable = false;
    }

    public void markAccessIssue() {
        accessIssue = true;
    }

    public boolean getAccessIssue() {
        return accessIssue;
    }

    public Type getTypeByInstanceCreation(ClassInstanceCreation instanceCreation) {
        return this.instanceCreationTypeMap.get(instanceCreation);
    }

    public boolean hasUnrefactorableNodePair() {
        return unrefactoredList.size() > 0;
    }

    public boolean comesFromSamePackage() {
        return packageDeclaration1.getName().getFullyQualifiedName().equals(
                packageDeclaration2.getName().getFullyQualifiedName());
    }

    public void addUnrefactoredNodePair(Expression node1, Expression node2, RFNodeDifference diff) {
        NodePair nodePair = new NodePair(node1, node2, diff);
        unrefactoredList.add(nodePair);
    }

    public void addTemplateArgumentPair(Expression arg1, Expression arg2) {
        this.templateArguments1.add((Expression) ASTNode.copySubtree(ast, arg1));
        this.templateArguments2.add((Expression) ASTNode.copySubtree(ast, arg2));
    }

    public boolean addVariableName(String name) {
        return templateVariableNames.add(name);
    }

    public void addStatement(Statement statement) {
        templateMethod.getBody().statements().add(statement);
    }

    public void addInstanceCreation(ClassInstanceCreation instanceCreation, Type type) {
        if (instanceCreation != null) {
            this.instanceCreationTypeMap.put(instanceCreation, type);
        }
    }

    public void addNullParameterObject() {
        if (!hasNullParameterObject) {
            VariableDeclarationFragment variableDeclarationFragment = ast.newVariableDeclarationFragment();
            variableDeclarationFragment.setName(ast.newSimpleName(NULL_PARAMETER_OBJECT));
            variableDeclarationFragment.setInitializer(ast.newNullLiteral());
            VariableDeclarationStatement variableDeclarationStatement = ast.newVariableDeclarationStatement(variableDeclarationFragment);
            variableDeclarationStatement.setType(ast.newSimpleType(ast.newSimpleName("Object")));
            templateMethod.getBody().statements().add(0, variableDeclarationStatement);
            hasNullParameterObject = true;
        }
    }

    public String resolveTypePair(TypePair typePair, boolean extendsCommonSuperClass) {
        if (!typeMap.containsKey(typePair)) {

            // get common super class
            ITypeBinding commonSuperClass = getLowestCommonSuperClass(typePair);
            ITypeBinding commonInterface = getLowestCommonInterface(typePair);

            //
            if (commonSuperClass == null && commonInterface == null
                    && (typePair.getType1().isParameterizedType() || typePair.getType2().isParameterizedType())) {
                markAsUnrefactorable();
                throw new IllegalStateException("no common super class for parameterized types: "
                        + typePair.getType1().getQualifiedName() + ", " + typePair.getType2().getQualifiedName());
            }

            // set type name
            String commonName = RenameUtil.constructCommonName(typePair.getType1().getName(),
                    typePair.getType2().getName(), true);
            if (commonInterface != null && commonInterface.isParameterizedType()) {
                commonName = commonInterface.getName();
            }
            commonName = RenameUtil.encodeParameterizedName(commonName);
            String typeName;
            if (!commonName.equals("")) {
                int nameCount = genericTypeNameMap.getOrDefault(commonName, 0);
                if (nameCount == 0) {
                    typeName = TYPE_NAME_PREFIX + commonName;
                } else {
                    typeName = TYPE_NAME_PREFIX + commonName + nameCount;
                }
                genericTypeNameMap.put(commonName, nameCount + 1);

            } else if (commonSuperClass != null) {

                String commonSuperClassName = commonSuperClass.getName();
                String[] nameComponents = RenameUtil.splitCamelCaseName(commonSuperClassName);
                StringBuilder newCommonSuperClassName = new StringBuilder();
                for (String component : nameComponents) {
                    if (!component.equals("Abstract")) {
                        newCommonSuperClassName.append(component);
                    }
                }

                int nameCount = genericTypeNameMap.getOrDefault(newCommonSuperClassName.toString(), 0);
                if (nameCount == 0) {
                    typeName = TYPE_NAME_PREFIX + newCommonSuperClassName.toString();
                } else {
                    typeName = TYPE_NAME_PREFIX + newCommonSuperClassName.toString() + nameCount;
                }
                genericTypeNameMap.put(newCommonSuperClassName.toString(), nameCount + 1);

            } else {
                typeName = TYPE_NAME_PREFIX + typeCount++;
            }
            typeMap.put(typePair, typeName);
            genericTypeMap.put(typeName, typePair);
            addGenericType(typeName);

            // add common generic type bound
            if (extendsCommonSuperClass) {
                if (commonSuperClass == null && commonInterface != null) {
                    commonSuperClass = commonInterface;
                }
                if (commonSuperClass != null) {
                    if (commonSuperClass.isParameterizedType()) {
                        List<TypeParameter> typeParameters = templateMethod.typeParameters();
                        for (TypeParameter typeParameter : typeParameters) {
                            if (typeParameter.getName().getIdentifier().equals(typeName)) {
                                typeParameter.typeBounds().add(ASTNodeUtil.typeFromBinding(ast, commonSuperClass));
                            }
                        }

                    } else {
                        addGenericTypeBound(typeName, commonSuperClass.getName());
                    }

                    String commonSuperClassPackageName = commonSuperClass.getPackage().getName();

                    // get common super class full name
                    String commonSuperClassFullName;
                    if (commonSuperClass.isParameterizedType()) {
                        commonSuperClassFullName = commonSuperClass.getErasure().getQualifiedName();
                    } else {
                        commonSuperClassFullName = commonSuperClass.getQualifiedName();
                    }

                    if (templateCU == null &&
                            !commonSuperClassPackageName.equals(packageDeclaration1.getName().getFullyQualifiedName())) {
                        templateCUImports.add(ASTNodeUtil.createPackageName(ast, commonSuperClassFullName));
                    }

                    if (templateCU != null &&
                            !commonSuperClassPackageName.equals(templateCU.getPackage().getName().getFullyQualifiedName())) {
                        addImportDeclaration(templateCU,
                                ASTNodeUtil.createPackageName(ast, commonSuperClassFullName), false);
                    }
                }
            }
        }
        return typeMap.get(typePair);
    }

    public static ITypeBinding getLowestCommonSuperClass(TypePair typePair) {
        ITypeBinding p1 = typePair.getType1();
        ITypeBinding p2 = typePair.getType2();

        while (p1 != null || p2 != null) {
            if (p1 != null && p2 != null && p1.getBinaryName().equals(p2.getBinaryName())) {
                return p1;

            } else {
                if (p1 == null) {
                    p1 = typePair.getType2();
                } else {
                    p1 = p1.getSuperclass();
                }

                if (p2 == null) {
                    p2 = typePair.getType1();
                } else {
                    p2 = p2.getSuperclass();
                }
            }
        }

        return null;
    }

    public static ITypeBinding getLowestCommonInterface(TypePair typePair) {
        ITypeBinding p1 = typePair.getType1();
        ITypeBinding p2 = typePair.getType2();

        Set<String> interfaceSet = new HashSet<>();
        Queue<ITypeBinding> queue1 = new ArrayDeque<>();
        Queue<ITypeBinding> queue2 = new ArrayDeque<>();

        if (p1.isInterface()) {
            queue1.offer(p1);
        } else {
            for (ITypeBinding iTypeBinding : p1.getInterfaces()) {
                queue1.offer(iTypeBinding);
            }
        }

        if (p2.isInterface()) {
            queue2.offer(p2);
        } else {
            for (ITypeBinding iTypeBinding : p2.getInterfaces()) {
                queue2.offer(iTypeBinding);
            }
        }

        while(!queue2.isEmpty()) {
            ITypeBinding cur = queue2.poll();
            interfaceSet.add(cur.getQualifiedName());
            for (ITypeBinding iTypeBinding : cur.getInterfaces()) {
                queue2.offer(iTypeBinding);
            }
        }

        while(!queue1.isEmpty()) {
            ITypeBinding cur = queue1.poll();
            if (interfaceSet.contains(cur.getQualifiedName())) {
                return cur;
            } else {
                for (ITypeBinding iTypeBinding : cur.getInterfaces()) {
                    queue1.offer(iTypeBinding);
                }
            }
        }
        return null;
    }

    public boolean containsVariableNamePair(String name1, String name2) {
        String resolvedName1 = nameMap1.getOrDefault(name1, "");
        String resolvedName2 = nameMap2.getOrDefault(name2, "");
        assert resolvedName1.equals(resolvedName2);
        return !resolvedName1.equals("");
    }

    public String resolveVariableName(String name1, String name2, String prefix) {
        String resolvedName1 = nameMap1.getOrDefault(name1, "");
        String resolvedName2 = nameMap2.getOrDefault(name2, "");
        assert resolvedName1.equals(resolvedName2);

        if (resolvedName1.equals("")) {
            String commonName = RenameUtil.renameVariable(name1, name2, variableCount++, prefix);
            while (!templateVariableNames.add(commonName)) {
                commonName = RenameUtil.renameVariable(name1, name2, variableCount++, prefix);
            }
            nameMap1.put(name1, commonName);
            nameMap2.put(name2, commonName);
            return commonName;

        } else {
            return resolvedName1;
        }
    }

    public String resolveGenericType(String genericType) {
        if (!clazzInstanceMap.containsKey(genericType)) {
            String clazzName = CLAZZ_NAME_PREFIX + genericType;
            //String clazzName = CLAZZ_NAME_PREFIX + genericType + clazzCount++;
            clazzInstanceMap.put(genericType, clazzName);
            recordLambdaForParameter(genericType);
            addClazzInParameter(genericType);
        }
        return clazzInstanceMap.get(genericType);
    }

    private void addTypeParameterAdapterImpl(Type type, List<Type> superInterfaceTypes) {
        Type interfaceType = superInterfaceTypes.get(0);
        ParameterizedType parameterizedType;
        if (interfaceType instanceof SimpleType) {
            parameterizedType =
                    ast.newParameterizedType((Type) ASTNode.copySubtree(ast, interfaceType));
            parameterizedType.typeArguments().add(type);
            superInterfaceTypes.remove(0);
            superInterfaceTypes.add(parameterizedType);

        } else {
            parameterizedType = (ParameterizedType) interfaceType;
            parameterizedType.typeArguments().add(type);
        }
    }

    private Type resolveAdapterActionArgumentType(Expression expr, ITypeBinding iTypeBinding) {

        // get exprType
        Type exprType;
        if (iTypeBinding != null && expr instanceof NullLiteral) {
            exprType = ASTNodeUtil.typeFromBinding(ast, iTypeBinding);
        } else {
            exprType = (Type) expr.getProperty(ASTNodeUtil.PROPERTY_TYPE_BINDING);
        }

        if (exprType != null) {
            if (typeMap.values().contains(exprType.toString())
                    && !adapterTypes.contains(exprType.toString())) {

                // add type parameter in adapter interface
                TypeParameter typeParameter = ast.newTypeParameter();
                typeParameter.setName(ast.newSimpleName(exprType.toString()));
                adapter.typeParameters().add(typeParameter);

                // add type parameter in adapter impl
                TypePair typePair = genericTypeMap.get(exprType.toString());
                if (typePair != null) {
                    Type type1 = ASTNodeUtil.typeFromBinding(ast, typePair.getType1());
                    Type type2 = ASTNodeUtil.typeFromBinding(ast, typePair.getType2());
                    addTypeParameterAdapterImpl(type1, adapterImpl1.superInterfaceTypes());
                    addTypeParameterAdapterImpl(type2, adapterImpl2.superInterfaceTypes());
                }

                // add adapter variable
                addAdapterVariableTypeParameter(exprType);

                adapterTypes.add(exprType.toString());
            }

        } else if (iTypeBinding != null) {
            exprType = ASTNodeUtil.typeFromBinding(ast, iTypeBinding);
        } else {
            exprType = ASTNodeUtil.typeFromBinding(ast, expr.resolveTypeBinding());
        }

        return exprType;
    }

    public void addImportDeclaration(CompilationUnit cu, Name name, boolean isStatic) {
        if (cu != null && name != null) {
            List<ImportDeclaration> importDeclarations = cu.imports();
            for (ImportDeclaration importDeclaration : importDeclarations) {
                if (importDeclaration.getName().getFullyQualifiedName().equals(name.getFullyQualifiedName())) {
                    return;
                }
            }
            if (name.isQualifiedName() && !isStatic) {
                QualifiedName qualifiedName = (QualifiedName) name;
                if (qualifiedName.getQualifier().getFullyQualifiedName().equals(cu.getPackage().getName().getFullyQualifiedName())) {
                    return;
                }
            }
            ImportDeclaration importDeclaration = ast.newImportDeclaration();
            importDeclaration.setName((Name) ASTNode.copySubtree(ast, name));
            importDeclaration.setStatic(isStatic);
            importDeclarations.add(importDeclaration);
        }
    }

    public void addImportInTemplateImportList(Name name) {
        if (name != null) {
            templateCUImports.add(name);
        }
    }

    public void addGenericTypeBound(String fromType, String toType) {
        List<TypeParameter> typeParameters = templateMethod.typeParameters();
        for (TypeParameter typeParameter : typeParameters) {
            if (typeParameter.getName().getIdentifier().equals(fromType)) {
                typeParameter.typeBounds().add(ast.newSimpleType(ast.newSimpleName(toType)));
                return;
            }
        }
    }

    private void addThrowsException() {
        if (templateMethod.thrownExceptionTypes().size() == 0) {
            templateMethod.thrownExceptionTypes().add(ast.newSimpleType(ast.newSimpleName(EXCEPTION_NAME)));
        }
    }

    private void addGenericType(String genericTypeName) {
        TypeParameter typeParameter = ast.newTypeParameter();
        typeParameter.setName(ast.newSimpleName(genericTypeName));
        templateMethod.typeParameters().add(typeParameter);
    }

    private void recordLambdaForParameter(String genericTypeName) {
        lambdaParameterType = computeLambdaVariableReturnType(genericTypeName);
    }

    private void actuallyAddLambdaParameter() {
        if (!lambdaParameterAdded && lambdaParameterType != null) {
            lambdaParameterAdded = true;
            addVariableParameter(lambdaParameterType, ast.newSimpleName(DEFAULT_LAMBDA_VARIABLE_NAME), 0);
        }
    }

    private void addClazzInParameter(String genericTypeName) {
        Type genericType = ast.newSimpleType(ast.newSimpleName(genericTypeName));
        Type classType = ast.newSimpleType(ast.newSimpleName(CLASS_NAME));
        ParameterizedType classTypeWithGenericType = ast.newParameterizedType(classType);
        classTypeWithGenericType.typeArguments().add(genericType);

        SimpleName clazzName = ast.newSimpleName(resolveGenericType(genericTypeName));
        addVariableParameter(classTypeWithGenericType, clazzName);

        TypePair typePair = genericTypeMap.get(genericTypeName);
        if (typePair != null) {
            TypeLiteral typeLiteral1 = ast.newTypeLiteral();
            Type type1;
            if (typePair.getType1().isParameterizedType()) {
                type1 = ASTNodeUtil.typeFromBinding(ast, typePair.getType1().getErasure());
            } else {
                type1 = ASTNodeUtil.typeFromBinding(ast, typePair.getType1());
            }
            typeLiteral1.setType(type1);
            templateArguments1.add(typeLiteral1);

            TypeLiteral typeLiteral2 = ast.newTypeLiteral();
            Type type2;
            if (typePair.getType2().isParameterizedType()) {
                type2 = ASTNodeUtil.typeFromBinding(ast, typePair.getType2().getErasure());
            } else {
                type2 = ASTNodeUtil.typeFromBinding(ast, typePair.getType2());
            }
            typeLiteral2.setType(type2);
            templateArguments2.add(typeLiteral2);
        }

        throwsAllException = true;
        //addThrowsException();
    }

    public String addVariableParameter(Type type) {
        int count = parameterMap.getOrDefault(type.toString(), 0) + 1;
        String variableParameter = RenameUtil.rename(type, count);
        while (!addVariableName(variableParameter)) {
            count++;
            variableParameter = RenameUtil.rename(type, count);
        }
        parameterMap.put(type.toString(), count);
        addVariableParameter(type, ast.newSimpleName(variableParameter));
        return variableParameter;
    }

    private void addVariableParameter(Type type, SimpleName name) {
        addVariableParameter(type, name, templateMethod.parameters().size());
    }

    private void addVariableParameter(Type type, SimpleName name, int index) {
        SingleVariableDeclaration variableParameter = ast.newSingleVariableDeclaration();
        variableParameter.setType(type);
        variableParameter.setName(name);
        templateMethod.parameters().add(index, variableParameter);
    }

    private void addAdapterVariableTypeParameter(Type type) {
        Type adapterType = adapterVariable.getType();
        if (adapterType.isSimpleType()) {
            ParameterizedType parameterizedType = ast.newParameterizedType((Type) ASTNode.copySubtree(ast, adapterType));
            parameterizedType.typeArguments().add(ASTNode.copySubtree(ast, type));
            adapterVariable.setType(parameterizedType);

        } else if (adapterType.isParameterizedType()) {
            ParameterizedType parameterizedType = (ParameterizedType) adapterType;
            parameterizedType.typeArguments().add(ASTNode.copySubtree(ast, type));
        } else {
            throw new IllegalStateException("unexpected adapter type");
        }
    }

    private void addMethodInAdapterInterface(SimpleName name, List<Type> argTypes, Type returnType,
                                             List<ITypeBinding> thrownExceptions, boolean isVarargs) {

        MethodDeclaration methodDeclaration = ast.newMethodDeclaration();

        // set return type
        methodDeclaration.setReturnType2((Type) ASTNode.copySubtree(ast, returnType));

        // add return type import declaration
        String returnTypeQualifiedName = (String) returnType.getProperty(ASTNodeUtil.PROPERTY_QUALIFIED_NAME);
        if (returnTypeQualifiedName != null) {
            addImportDeclaration(adapterInterfaceCU, ASTNodeUtil.createPackageName(ast, returnTypeQualifiedName), false);
        }

        // set interface action name
        methodDeclaration.setName((SimpleName) ASTNode.copySubtree(ast, name));

        Map<String, Integer> argMap = new HashMap<>();
        for (int i = 0; i < argTypes.size(); i++) {

            Type argType = argTypes.get(i);
            if (argType == null) continue;

            SingleVariableDeclaration arg = ast.newSingleVariableDeclaration();

            // set arg type
            if (i == argTypes.size() - 1 && isVarargs) {
                ArrayType arrayType = (ArrayType) argType;
                arg.setType((Type) ASTNode.copySubtree(ast, arrayType.getElementType()));
                arg.setVarargs(true);
            } else {
                arg.setType((Type) ASTNode.copySubtree(ast, argType));
            }

            // add import declaration
            String qualifiedName = (String) argType.getProperty(ASTNodeUtil.PROPERTY_QUALIFIED_NAME);
            if (qualifiedName != null && templateCU != null) {
                addImportDeclaration(adapterInterfaceCU, ASTNodeUtil.createPackageName(ast, qualifiedName), false);
            }

            // set arg name
            String argName = RenameUtil.rename(argType, 0);
            int argCount = argMap.getOrDefault(argName, 1);
            argMap.put(argName, argCount + 1);
            arg.setName(ast.newSimpleName(RenameUtil.rename(argType, argCount)));

            // add parameter
            methodDeclaration.parameters().add(arg);
        }

        // add throw exception
        for (ITypeBinding throwsException : thrownExceptions) {
            addThrownExceptionTypesInMethodDeclaration(methodDeclaration, throwsException);
            addImportDeclaration(templateCU,
                    ASTNodeUtil.createPackageName(ast, throwsException.getBinaryName()), false);
        }

        adapter.bodyDeclarations().add(methodDeclaration);
    }

    private void addThrownExceptionTypesInMethodDeclaration(MethodDeclaration methodDeclaration, ITypeBinding throwsException) {
        List<Type> thrownExceptionTypes = methodDeclaration.thrownExceptionTypes();
        for (Type exceptionType : thrownExceptionTypes) {
            ITypeBinding cur = throwsException;
            while (cur != null) {
                if (cur.getName().equals(exceptionType.toString())) {
                    return;
                }
                cur = cur.getSuperclass();
            }
        }
        thrownExceptionTypes.add(ASTNodeUtil.typeFromBinding(ast, throwsException));
    }

    private MethodDeclaration addMethodInAdapterImpl(SimpleName actionName, List<Type> argTypes,
                                                     MethodInvocationPair methodInvocationPair,
                                                     ITypeBinding returnTypeBinding, Pair pair) {

        Expression expr;
        SimpleName name;
        List<Expression> arguments;
        IMethodBinding iMethodBinding;
        MethodInvocation curMethodInvocation;
        ITypeBinding[] parameterTypeBinding;

        // get cu
        CompilationUnit cu;
        if (pair == Pair.member1) {
            cu = adapterImplCU1;
        } else {
            cu = adapterImplCU2;
        }

        switch (pair) {
            case member1:
                expr = methodInvocationPair.getExpr1();
                name = methodInvocationPair.getName1();
                arguments = methodInvocationPair.getArgument1();
                iMethodBinding = methodInvocationPair.getMethod1().resolveMethodBinding();
                curMethodInvocation = methodInvocationPair.getMethod1();
                parameterTypeBinding = methodInvocationPair.getExtendArgTypeBinding1();
                break;
            case member2:
            default:
                expr = methodInvocationPair.getExpr2();
                name = methodInvocationPair.getName2();
                arguments = methodInvocationPair.getArgument2();
                iMethodBinding = methodInvocationPair.getMethod2().resolveMethodBinding();
                curMethodInvocation = methodInvocationPair.getMethod2();
                parameterTypeBinding = methodInvocationPair.getExtendArgTypeBinding2();
        }

        MethodDeclaration method = ast.newMethodDeclaration();
        method.setBody(ast.newBlock());

        // set return type
        Type returnType;
        if (returnTypeBinding != null) {
            if (!returnTypeBinding.isPrimitive() && !returnTypeBinding.getBinaryName().startsWith(JAVA_LANG_CLASS)) {
                addImportDeclaration(cu, ASTNodeUtil.createPackageName(ast, returnTypeBinding.getBinaryName()), false);
            }
            returnType = ASTNodeUtil.typeFromBinding(ast, returnTypeBinding);
        } else {
            returnType = ast.newSimpleType(ast.newSimpleName(OBJECT_NAME));
            returnType.setProperty(ASTNodeUtil.PROPERTY_QUALIFIED_NAME, JAVA_OBJECT_FULL_NAME);
        }
        method.setReturnType2((Type) ASTNode.copySubtree(ast, returnType));

        // set interface action name
        method.setName((SimpleName) ASTNode.copySubtree(ast, actionName));

        // set modifier
        method.modifiers().add(ast.newModifier(Modifier.ModifierKeyword.PUBLIC_KEYWORD));

        // create method invocation
        MethodInvocation methodInvocation = ast.newMethodInvocation();
        methodInvocation.setName((SimpleName) ASTNode.copySubtree(ast, name));
        Map<String, Integer> argMap = new HashMap<>();

        for (int i = 0; i < argTypes.size(); i++) {

            // set method variable declaration type
            SingleVariableDeclaration variableDeclaration = ast.newSingleVariableDeclaration();
            Type curType = argTypes.get(i);
            if (i ==  argTypes.size() - 1 && iMethodBinding.isVarargs()) {
                curType = ((ArrayType) curType).getElementType();
                variableDeclaration.setVarargs(true);
            }
            if (curType == null) {
                if (expr instanceof SimpleName && ((SimpleName) expr).resolveBinding().getKind() == 2) {
                    // Type static call
                    methodInvocation.setExpression((Expression) ASTNode.copySubtree(ast, expr));
                }
                continue;
            }
            String argTypeName = curType.toString();

            Type argType;
            if (genericTypeMap.containsKey(argTypeName)) {
                TypePair typePair = genericTypeMap.get(argTypeName);

                //Type argType;
                if (pair == Pair.member1) {
                    argType = ASTNodeUtil.typeFromBinding(ast, typePair.getType1());
                } else {
                    argType = ASTNodeUtil.typeFromBinding(ast, typePair.getType2());
                }

            } else {
                argType = ASTNodeUtil.copyTypeWithProperties(ast, curType);
            }
            variableDeclaration.setType(argType);

            // add import declaration
            addImportDeclaration(cu,
                    ASTNodeUtil.createPackageName(ast, (String) argType.getProperty(ASTNodeUtil.PROPERTY_QUALIFIED_NAME)), false);

            // get current expression
            Expression curExpr;
            if (i == 0) {
                curExpr = expr;
            } else {
                curExpr = arguments.get(i - 1);
            }

            // set method expr variable name
            SimpleName curName;
            if (curExpr instanceof SimpleName) {
                SimpleName curExprSimpleName = (SimpleName) curExpr;
                if (argMap.containsKey(curExprSimpleName.getIdentifier())) {
                    int argCount = argMap.get(curExprSimpleName.getIdentifier()) + 1;
                    curName = ast.newSimpleName(curExprSimpleName.getIdentifier() + argCount);
                    argMap.put(curExprSimpleName.getIdentifier(), argCount);

                } else {
                    curName = curExprSimpleName;
                    argMap.put(curExprSimpleName.getIdentifier(), 0);
                }

            } else {
                String argName = RenameUtil.rename(curType, 0);
                int argCount = argMap.getOrDefault(argName, 1);
                curName = ast.newSimpleName(RenameUtil.rename(curType, argCount));
                argMap.put(argName, argCount + 1);
            }
            variableDeclaration.setName((SimpleName) ASTNode.copySubtree(ast, curName));
            method.parameters().add(variableDeclaration);

            if (i == 0) {
                methodInvocation.setExpression((Expression) ASTNode.copySubtree(ast, curName));

            } else {
                // check if argType matches iMethodBinding
                ITypeBinding parameterType = parameterTypeBinding[i - 1];
                String qualifiedName = (String) argType.getProperty(ASTNodeUtil.PROPERTY_QUALIFIED_NAME);

                if (qualifiedName != null && !parameterType.getBinaryName().equals(qualifiedName)) {
                    // check if argType is subType of parameterType
                    if (ASTNodeUtil.hasAncestor(parameterType, qualifiedName)) {
                        if (!parameterType.isPrimitive()) {
                            addImportDeclaration(cu,
                                    ASTNodeUtil.createPackageName(ast, parameterType.getBinaryName()), false);
                        }
                        CastExpression castExpression = ast.newCastExpression();
                        castExpression.setExpression((Expression) ASTNode.copySubtree(ast, curName));
                        castExpression.setType(ASTNodeUtil.typeFromBinding(ast, parameterType));
                        methodInvocation.arguments().add(castExpression);
                    } else {
                        methodInvocation.arguments().add(ASTNode.copySubtree(ast, curName));
                    }

                } else {
                    methodInvocation.arguments().add(ASTNode.copySubtree(ast, curName));
                }
            }

        }

        // add statement to method
        Statement statement;
        if (returnType.isPrimitiveType() && ((PrimitiveType) returnType).getPrimitiveTypeCode() == PrimitiveType.VOID) {
            statement = ast.newExpressionStatement(methodInvocation);

        } else {

            ReturnStatement returnStatement = ast.newReturnStatement();

            // check if return type is compatible
            CastExpression castExpression = ast.newCastExpression();
            castExpression.setType(ASTNodeUtil.copyTypeWithProperties(ast, returnType));

            boolean isSameType = true;
            if (returnTypeBinding != null && curMethodInvocation.resolveTypeBinding() != null) {

                if (returnTypeBinding.isParameterizedType() && curMethodInvocation.resolveTypeBinding().isParameterizedType()) {
                    ITypeBinding[] typeArgs1 = returnTypeBinding.getTypeArguments();
                    ITypeBinding[] typeArgs2 = curMethodInvocation.resolveTypeBinding().getTypeArguments();
                    if (typeArgs1.length != typeArgs2.length) {
                        isSameType = false;
                    } else {
                        for (int i = 0; i < typeArgs1.length; i++) {
                            ITypeBinding typeArg1 = typeArgs1[i];
                            ITypeBinding typeArg2 = typeArgs2[i];
                            if (typeArg1.isCapture()) {
                                typeArg1 = typeArg1.getWildcard();
                            }
                            if (typeArg2.isCapture()) {
                                typeArg2 = typeArg2.getWildcard();
                            }
                            if (!typeArg1.getQualifiedName().equals(typeArg2.getQualifiedName())) {
                                isSameType = false;
                                break;
                            }
                        }
                    }

                } else {
                    if (!returnTypeBinding.getQualifiedName().equals(curMethodInvocation.resolveTypeBinding().getQualifiedName())
                            && !curMethodInvocation.resolveTypeBinding().isAssignmentCompatible(returnTypeBinding)) {
                        isSameType = false;
                    }
                }
            }

            if (isSameType) {
                returnStatement.setExpression(methodInvocation);
            } else {
                castExpression.setExpression(methodInvocation);
                returnStatement.setExpression(castExpression);
            }
            statement = returnStatement;
        }
        method.getBody().statements().add(statement);

        // add throws exception
        for (ITypeBinding thrownException : iMethodBinding.getExceptionTypes()) {
            addThrownExceptionTypesInMethodDeclaration(method, thrownException);
        }

        return method;
    }

    private MethodDeclaration addMethodInAdapterImpl(SimpleName actionName, List<Type> argTypes, Type returnType, Pair pair) {

        MethodDeclaration method = ast.newMethodDeclaration();
        method.setBody(ast.newBlock());

        // set return type
        method.setReturnType2((Type) ASTNode.copySubtree(ast, returnType));

        // set interface action name
        method.setName((SimpleName) ASTNode.copySubtree(ast, actionName));

        // set modifier
        method.modifiers().add(ast.newModifier(Modifier.ModifierKeyword.PUBLIC_KEYWORD));

        assert argTypes.size() == 2;

        // get cu
        CompilationUnit cu;
        if (pair == Pair.member1) {
            cu = adapterImplCU1;
        } else {
            cu = adapterImplCU2;
        }

        // set method parameters
        Map<String, Integer> argMap = new HashMap<>();
        for (Type argType : argTypes) {
            SingleVariableDeclaration variableDeclaration = ast.newSingleVariableDeclaration();
            variableDeclaration.setType((Type) ASTNode.copySubtree(ast, argType));

            // add import declaration
            if (templateCU != null) {
                addImportDeclaration(cu,
                        ASTNodeUtil.createPackageName(ast, (String) argType.getProperty(ASTNodeUtil.PROPERTY_QUALIFIED_NAME)), false);
            }

            int argCount = argMap.getOrDefault(argType.toString(), 1);
            variableDeclaration.setName(ast.newSimpleName(RenameUtil.rename(argType, argCount)));
            argMap.put(argType.toString(), argCount + 1);
            method.parameters().add(variableDeclaration);
        }

        // get return variable name
        SimpleName returnVariableName;
        if (pair == Pair.member1) {
            returnVariableName = ((SingleVariableDeclaration) method.parameters().get(0)).getName();
        } else {
            returnVariableName = ((SingleVariableDeclaration) method.parameters().get(1)).getName();
        }

        // add return statement to method
        ReturnStatement returnStatement = ast.newReturnStatement();
        returnStatement.setExpression((Expression) ASTNode.copySubtree(ast, returnVariableName));
        method.getBody().statements().add(returnStatement);

        return method;
    }

    private void addAdapterActionImpl(SimpleName actionName, List<Type> argTypes, MethodInvocationPair pair,
                                      TypePair returnTypePair) {
        ITypeBinding returnTypeBinding1;
        ITypeBinding returnTypeBinding2;
        TypePair boundPair = ASTNodeUtil.getJavaLangClassCaptureBounds(returnTypePair);
        if (containsTypePair(returnTypePair) || (boundPair != null && containsTypePair(boundPair))) {
            returnTypeBinding1 = returnTypePair.getType1();
            returnTypeBinding2 = returnTypePair.getType2();
        } else {
            ITypeBinding returnTypeBinding = ASTNodeUtil.getAssignmentCompatibleTypeBinding(returnTypePair);
            returnTypeBinding1 = returnTypeBinding;
            returnTypeBinding2 = returnTypeBinding;
        }

        MethodDeclaration method1 = addMethodInAdapterImpl(actionName, argTypes, pair, returnTypeBinding1, Pair.member1);
        MethodDeclaration method2 = addMethodInAdapterImpl(actionName, argTypes, pair, returnTypeBinding2, Pair.member2);
        adapterImpl1.bodyDeclarations().add(method1);
        adapterImpl2.bodyDeclarations().add(method2);
    }

    private void addAdapterActionImpl(SimpleName actionName, List<Type> argTypes, Type returnType) {
        MethodDeclaration method1 = addMethodInAdapterImpl(actionName, argTypes, returnType, Pair.member1);
        MethodDeclaration method2 = addMethodInAdapterImpl(actionName, argTypes, returnType, Pair.member2);
        adapterImpl1.bodyDeclarations().add(method1);
        adapterImpl2.bodyDeclarations().add(method2);
    }

    private void addAdapterVariableParameter() {
        if (!hasAdapterVariable) {
            templateMethod.parameters().add(0, adapterVariable);
            hasAdapterVariable = true;
        }
    }

    private void addLambda() {
        if (ENABLE_LAMBDA) {
            hasLambda = true;
        }
    }

    private boolean containsTypeVariable(ITypeBinding iTypeBinding) {
        if (iTypeBinding == null) {
            return false;
        }

        if (iTypeBinding.isArray()) {
            return containsTypeVariable(iTypeBinding.getElementType());

        } else if (iTypeBinding.isParameterizedType()) {
            for (ITypeBinding typeArg: iTypeBinding.getTypeArguments()) {
                if (containsTypeVariable(typeArg)) {
                    return true;
                }
            }
            return containsTypeVariable(iTypeBinding.getErasure());

        } else if (iTypeBinding.isWildcardType()) {
            return containsTypeVariable(iTypeBinding.getBound());
        }

        return iTypeBinding.isTypeVariable();
    }

    private Type computeLambdaVariableReturnType(String rt) {
        String olTypeName = String.format(DEFAULT_LAMBDA_FUNCTION_CLASS_TEMPLATE, "Long");
        Type openLambdaType = ast.newSimpleType(ast.newName(olTypeName));
        ParameterizedType lType = ast.newParameterizedType(openLambdaType);
        Type t = ast.newSimpleType(ast.newSimpleName(rt));
        lType.typeArguments().add(t);
        return lType;
    }

    @SuppressWarnings("unchecked")
    public MethodInvocation createLambdaActionMethodInvocation(Expression expr, List<Expression> arguments,
            MethodInvocationPair pair, TypePair returnTypePair) {
        // create new method invocation
        MethodInvocation newMethodInvocation = ast.newMethodInvocation();
        newMethodInvocation.setExpression(ast.newSimpleName(lambdaVariable.getName().getIdentifier()));

        List<Expression> newArgs = newMethodInvocation.arguments();
        List<Type> argTypes = new ArrayList<>();

        // XXX PL need to handle lambdas that return other than int, long, double
        Type returnType = computeReturnTypeFromPair(returnTypePair);
        boolean foundWinningType = false;

        if (returnType.isPrimitiveType()) {
            PrimitiveType pt = (PrimitiveType) returnType;
            if (pt.getPrimitiveTypeCode() == PrimitiveType.INT) {
                newMethodInvocation.setName(ast.newSimpleName("applyAsInt"));
                lambdaVariable.setType(computeLambdaVariableReturnType("Int"));
                foundWinningType = true;
            } else if (pt.getPrimitiveTypeCode() == PrimitiveType.LONG) {
                newMethodInvocation.setName(ast.newSimpleName("applyAsLong"));
                lambdaVariable.setType(computeLambdaVariableReturnType("Long"));
                foundWinningType = true;
            } else if (pt.getPrimitiveTypeCode() == PrimitiveType.DOUBLE) {
                newMethodInvocation.setName(ast.newSimpleName("applyAsDouble"));
                lambdaVariable.setType(computeLambdaVariableReturnType("Double"));
                foundWinningType = true;
            }
        }
        if (!foundWinningType) {
            log.info("could not do a lambda with return type "+returnType);
        }

        if (expr != null) {
            newArgs.add((Expression) ASTNode.copySubtree(ast, expr));
            argTypes.add(resolveAdapterActionArgumentType(expr, null));
        }

        if (foundWinningType) {
            addLambda();
            actuallyAddLambdaParameter();
            return newMethodInvocation;
        } else
            return null;
    }

    public MethodInvocation createAdapterActionMethod(Expression expr, List<Expression> arguments,
                                                      MethodInvocationPair pair, TypePair returnTypePair) {

        boolean ENABLE_THIS_LAMBDA = true;
        // XXX PL need to figure out when lambdas should be created vs adapters...
        if (ENABLE_LAMBDA && ENABLE_THIS_LAMBDA) {
            // stash away (in global state) the body --- which expr we're supposed to put in place of method1/method2
            method1.setProperty(ASTNodeUtil.PROPERTY_LAMBDA_PARAM, pair.getMethod1());
            method2.setProperty(ASTNodeUtil.PROPERTY_LAMBDA_PARAM, pair.getMethod2());

            MethodInvocation mi = createLambdaActionMethodInvocation(expr, arguments, pair, returnTypePair);
            if (mi != null)
                return mi;
        }

        addAdapterVariableParameter();

        // create new method invocation
        MethodInvocation newMethod = ast.newMethodInvocation();
        newMethod.setExpression(ast.newSimpleName(adapterVariable.getName().getIdentifier()));

        List<Expression> newArgs = newMethod.arguments();
        List<Type> argTypes = new ArrayList<>();

        // check expr
        if (expr != null) {
            Type exprType = (Type) expr.getProperty(ASTNodeUtil.PROPERTY_TYPE_BINDING);
            String methodName1 = pair.getMethod1().getName().getIdentifier();
            String methodName2 = pair.getMethod2().getName().getIdentifier();
            if (exprType != null && JAVA_OBJECT_FULL_NAME.equals(exprType.getProperty(ASTNodeUtil.PROPERTY_QUALIFIED_NAME))
                    && (!javaObjectMethodNames.contains(methodName1) || !javaObjectMethodNames.contains(methodName2))) {
                markAsUnrefactorable();
                log.info("cannot refactor incompatible method calls in java.lang.Object: " + methodName1 + ", " + methodName2);
            }
        }

        // copy and resolve method expr
        if (expr != null) {
            if (expr instanceof SimpleName && ((SimpleName) expr).resolveBinding().getKind() == 2) {
                argTypes.add(null);

            } else {
                newArgs.add((Expression) ASTNode.copySubtree(ast, expr));
                argTypes.add(resolveAdapterActionArgumentType(expr, null));
            }
        } else {
            // mark adapter impl as inner class to make use of class private methods
            //innerImplClass = true;
            argTypes.add(null);
        }


        // copy and resolve arguments
        for (int i = 0; i < arguments.size(); i++) {
            Expression argument = arguments.get(i);
            newArgs.add((Expression) ASTNode.copySubtree(ast, argument));
        }

        ITypeBinding[] iTypeBindings1 = pair.getMethod1().resolveMethodBinding().getParameterTypes();
        ITypeBinding[] iTypeBindings2 = pair.getMethod2().resolveMethodBinding().getParameterTypes();
        boolean isVarargs = pair.getMethod1().resolveMethodBinding().isVarargs();
        for (int i = 0; i < iTypeBindings1.length; i++) {
            if (containsTypeVariable(iTypeBindings1[i]) || containsTypeVariable(iTypeBindings2[i])) {
                log.info("unable to refactor type pair ["
                        + iTypeBindings1[i].getQualifiedName() + ", " + iTypeBindings2[i].getQualifiedName()
                        + "] with type variable when creating adapter action method.");
                markAsUnrefactorable();
            }
            if (iTypeBindings1[i].getQualifiedName().equals(iTypeBindings2[i].getQualifiedName())) {
                argTypes.add(ASTNodeUtil.typeFromBinding(ast, iTypeBindings1[i]));
            } else {
                argTypes.add(resolveAdapterActionArgumentType(arguments.get(i), iTypeBindings1[i]));
            }
        }

        // add method in adapter interface
        if (methodInvocationMap.containsKey(pair)) {
            String methodName = methodInvocationMap.get(pair);
            newMethod.setName(ast.newSimpleName(methodName));
            for (MethodDeclaration methodDeclaration : adapter.getMethods()) {
                if (methodDeclaration.getName().getIdentifier().equals(methodName)) {
                    Type returnType = methodDeclaration.getReturnType2();
                    newMethod.setProperty(ASTNodeUtil.PROPERTY_TYPE_BINDING, ASTNodeUtil.copyTypeWithProperties(ast, returnType));
                }
            }

        } else {

            // set adapter action name
            String newActionName;
            if (pair.getName1() != null && pair.getName2() != null) {
                String commonName = RenameUtil.constructCommonName(pair.getName1().getIdentifier(),
                        pair.getName2().getIdentifier(), false);
                if (!commonName.equals("")) {

                    if (keywords.contains(commonName)) {
                        commonName += KEYWORD_SUFFIX;
                    }

                    int nameCount = adapterActionNameMap.getOrDefault(commonName, 0);
                    if (nameCount == 0) {
                        newActionName = commonName;
                    } else {
                        newActionName = commonName + nameCount;
                    }
                    adapterActionNameMap.put(commonName, nameCount + 1);

                } else {
                    newActionName = DEFAULT_ADAPTER_METHOD_NAME + actionCount++;
                }
            } else {
                newActionName = DEFAULT_ADAPTER_METHOD_NAME + actionCount++;
            }
            newMethod.setName(ast.newSimpleName(newActionName));

            // add method in adapter interface
            Type returnType = computeReturnTypeFromPair(returnTypePair);

            // check if it needs to throw exception
            IMethodBinding methodBinding1 = pair.getMethod1().resolveMethodBinding();
            IMethodBinding methodBinding2 = pair.getMethod2().resolveMethodBinding();
            List<ITypeBinding> thrownExceptions = new ArrayList<>();
            for (ITypeBinding iTypeBinding : methodBinding1.getExceptionTypes()) {
                thrownExceptions.add(iTypeBinding);
            }
            for (ITypeBinding iTypeBinding : methodBinding2.getExceptionTypes()) {
                thrownExceptions.add(iTypeBinding);
            }

            addMethodInAdapterInterface(newMethod.getName(), argTypes, returnType, thrownExceptions, isVarargs);
            newMethod.setProperty(ASTNodeUtil.PROPERTY_TYPE_BINDING, ASTNodeUtil.copyTypeWithProperties(ast, returnType));
            methodInvocationMap.put(pair, newActionName);

            // create adapter action impl
            addAdapterActionImpl(newMethod.getName(), argTypes, pair, returnTypePair);

            // add throws exception in template method
            for (ITypeBinding thrownException : thrownExceptions) {
                addThrownExceptionTypesInMethodDeclaration(templateMethod, thrownException);
            }
        }

        return newMethod;
    }

    private Type computeReturnTypeFromPair(TypePair returnTypePair) {
        ITypeBinding returnTypeBinding = ASTNodeUtil.getAssignmentCompatibleTypeBinding(returnTypePair);
        TypePair boundPair = ASTNodeUtil.getJavaLangClassCaptureBounds(returnTypePair);
        Type returnType;

        if (containsTypePair(returnTypePair)) {

            // add type parameter in adapter interface
            String genericName = resolveTypePair(returnTypePair, false);
            returnType = ast.newSimpleType(ast.newSimpleName(genericName));
            if (!adapterTypes.contains(genericName)) {
                TypeParameter typeParameter = ast.newTypeParameter();
                typeParameter.setName(ast.newSimpleName(genericName));
                adapter.typeParameters().add(typeParameter);

                // add type parameter in adapter impl
                Type type1 = ASTNodeUtil.typeFromBinding(ast, returnTypePair.getType1());
                Type type2 = ASTNodeUtil.typeFromBinding(ast, returnTypePair.getType2());
                addTypeParameterAdapterImpl(type1, adapterImpl1.superInterfaceTypes());
                addTypeParameterAdapterImpl(type2, adapterImpl2.superInterfaceTypes());

                // add adapter variable
                addAdapterVariableTypeParameter(returnType);

                adapterTypes.add(genericName);
            }

        } else if (boundPair != null && containsTypePair(boundPair)) {
            String genericName = resolveTypePair(boundPair, false);
            ParameterizedType parameterizedType =
                    ast.newParameterizedType(ASTNodeUtil.typeFromBinding(ast, returnTypePair.getType1().getErasure()));
            WildcardType capType = ast.newWildcardType();
            capType.setBound(ast.newSimpleType(ast.newSimpleName(genericName)), true);
            parameterizedType.typeArguments().add(capType);
            returnType = parameterizedType;

        } else {
            if (returnTypeBinding == null) {
                returnType = ast.newSimpleType(ast.newSimpleName(OBJECT_NAME));
                returnType.setProperty(ASTNodeUtil.PROPERTY_QUALIFIED_NAME, JAVA_OBJECT_FULL_NAME);

            } else {
                if (!returnTypeBinding.isPrimitive() && !returnTypeBinding.getBinaryName().startsWith(JAVA_LANG_CLASS)) {
                    addImportDeclaration(templateCU, ASTNodeUtil.createPackageName(ast, returnTypeBinding.getBinaryName()), false);
                }
                returnType = ASTNodeUtil.typeFromBinding(ast, returnTypeBinding);
            }
        }
        return returnType;
    }

    public MethodInvocation createAdapterActionMethod(Expression e1, Expression e2, Type returnType) {

        addAdapterVariableParameter();

        MethodInvocation newMethod = ast.newMethodInvocation();
        newMethod.setExpression(ast.newSimpleName(adapterVariable.getName().getIdentifier()));

        String newActionName = DEFAULT_ADAPTER_METHOD_NAME + actionCount++;
        newMethod.setName(ast.newSimpleName(newActionName));

        List<Expression> newArgs = newMethod.arguments();
        List<Type> argTypes = new ArrayList<>();
        newArgs.add((Expression) ASTNode.copySubtree(ast, e1));
        newArgs.add((Expression) ASTNode.copySubtree(ast, e2));
        argTypes.add(resolveAdapterActionArgumentType(e1, null));
        argTypes.add(resolveAdapterActionArgumentType(e2, null));

        // add method in adapter interface
        addMethodInAdapterInterface(newMethod.getName(), argTypes, returnType, new ArrayList<>(), false);

        // create adapter action impl
        addAdapterActionImpl(newMethod.getName(), argTypes, returnType);

        return newMethod;
    }

    /*
    public MethodInvocation createAdapterActionMethod(Type returnType) {

        addAdapterVariableParameter();

        MethodInvocation newMethod = ast.newMethodInvocation();
        newMethod.setExpression(ast.newSimpleName(adapterVariable.getName().getIdentifier()));

        String newActionName = DEFAULT_ADAPTER_METHOD_NAME + actionCount++;
        newMethod.setName(ast.newSimpleName(newActionName));
        addMethodInAdapterInterface(newMethod.getName(), new ArrayList<>(), returnType);

        return newMethod;
    }
    */

    private boolean methodContainsThrowable() {
        List<Type>  thrownExceptionTypes1 = method1.thrownExceptionTypes();
        List<Type>  thrownExceptionTypes2 = method2.thrownExceptionTypes();
        for (Type thrownExceptionType : thrownExceptionTypes1) {
            if (thrownExceptionType.toString().equals(THROWABLE_NAME)) {
                return true;
            }
        }

        for (Type thrownExceptionType : thrownExceptionTypes2) {
            if (thrownExceptionType.toString().equals(THROWABLE_NAME)) {
                return true;
            }
        }

        return false;
    }

    public void modifyTestMethods() {
        if (throwsAllException || templateMethod.thrownExceptionTypes().size() > 1) {
            throwsAllException = true;
            templateMethod.thrownExceptionTypes().clear();
            if (methodContainsThrowable()) {
                templateMethod.thrownExceptionTypes().add(ast.newSimpleType(ast.newSimpleName(THROWABLE_NAME)));
            } else {
                templateMethod.thrownExceptionTypes().add(ast.newSimpleType(ast.newSimpleName(EXCEPTION_NAME)));
            }
        }
        modifyMethod(method1, adapterImpl1, templateArguments1, Pair.member1);
        modifyMethod(method2, adapterImpl2, templateArguments2, Pair.member2);
    }

    /** swaps out the implementation of method body with a call to templateMethod */
    private void modifyMethod(MethodDeclaration method, TypeDeclaration adapterImpl, List<Expression> arguments, Pair pair) {
        // create new method invocation
        MethodInvocation methodInvocation = ast.newMethodInvocation();
        methodInvocation.setName((SimpleName) ASTNode.copySubtree(ast, templateMethod.getName()));
        if (templateClass != null) {
            methodInvocation.setExpression(ast.newSimpleName(templateClass.getName().getIdentifier()));
        } else {
            methodInvocation.setExpression(ast.newThisExpression());
        }

        // add method arguments
        List<Expression> args = methodInvocation.arguments();

        // add type parameter in template method invocation
        if (!hasAdapterVariable && arguments.isEmpty()) {
            List<TypeParameter> typeParameters = templateMethod.typeParameters();
            for (TypeParameter typeParameter : typeParameters) {
                TypePair typePair = genericTypeMap.get(typeParameter.getName().getIdentifier());
                if (typePair != null) {
                    Type type;
                    if (pair == Pair.member1) {
                        type = ASTNodeUtil.typeFromBinding(ast, typePair.getType1());
                    } else {
                        type = ASTNodeUtil.typeFromBinding(ast, typePair.getType2());
                    }
                    methodInvocation.typeArguments().add(type);
                }
            }
        }

        if (ENABLE_LAMBDA && hasLambda) {
            // currently support single-arg lambdas that return int, long, or double
            LambdaExpression lambdaExpression = ast.newLambdaExpression();
            // XXX PL args and body
            // was: target.get()
            // want: (AtomicInteger target) -> target.get();

            VariableDeclarationFragment arg = ast.newVariableDeclarationFragment();
            arg.setName(ast.newSimpleName(DEFAULT_LAMBDA_VARIABLE_NAME));
            lambdaExpression.parameters().add(arg);

            MethodInvocation body = ast.newMethodInvocation();
            if (method.properties().containsKey(ASTNodeUtil.PROPERTY_LAMBDA_PARAM)) {
                body = (MethodInvocation) ASTNode.copySubtree(ast,
                        (ASTNode)method.getProperty(ASTNodeUtil.PROPERTY_LAMBDA_PARAM));
                body.setExpression(ast.newSimpleName(DEFAULT_LAMBDA_VARIABLE_NAME));
            }

            lambdaExpression.setBody(body);

            args.add(lambdaExpression);
        }

        if (hasAdapterVariable) {
            ClassInstanceCreation classInstanceCreation = ast.newClassInstanceCreation();
            SimpleName simpleName = (SimpleName) ASTNode.copySubtree(ast, adapterImpl.getName());
            classInstanceCreation.setType(ast.newSimpleType(simpleName));
            args.add(classInstanceCreation);
        }

        for (Expression arg : arguments) {
            args.add((Expression) ASTNode.copySubtree(ast, arg));
        }

        // add method invocation to method body
        Block body = ast.newBlock();
        body.statements().add(ast.newExpressionStatement(methodInvocation));
        method.setBody(body);

        // rename methods
        //method.getName().setIdentifier(method.getName().getIdentifier() + "RF");

        // remove javadoc
        //method.setJavadoc(null);

        // add Exception handling if necessary
        if (throwsAllException) {
            method.thrownExceptionTypes().clear();
            if (methodContainsThrowable()) {
                method.thrownExceptionTypes().add(ast.newSimpleType(ast.newSimpleName(THROWABLE_NAME)));
            } else {
                method.thrownExceptionTypes().add(ast.newSimpleType(ast.newSimpleName(EXCEPTION_NAME)));
            }
        }
    }

    public String getTemplatName() {
        return templateMethod.getName().getFullyQualifiedName();
    }

    public void updateSourceFiles() throws Exception {

        log.info("refactoring method pair: "
                + packageDeclaration1.getName().getFullyQualifiedName() + "."
                + compilationUnit1.getJavaElement().getElementName()
                + "#" + method1.getName().getIdentifier() + " <---> "
                + packageDeclaration2.getName().getFullyQualifiedName() + "."
                + compilationUnit2.getJavaElement().getElementName()
                + "#" + method2.getName().getIdentifier());

        if (templateCU == null) {
            // duplicate methods are in the same file
            saveMethod(iCU1, templateMethod, templateCUImports, false);
            if (hasAdapterVariable) {
                insertInnerClass(iCU1, adapterInterfaceCU, true);
                insertInnerClass(iCU1, adapterImplCU1, true);
                insertInnerClass(iCU1, adapterImplCU2, true);
            }
            cleanImportDeclarations(iCU1);

        } else if (packageDeclaration1.getName().getFullyQualifiedName().equals(packageDeclaration2.getName().getFullyQualifiedName())) {
            // duplicate methods are in the same package but different files
            ICompilationUnit iCompilationUnit = saveCU(packageFragmentRoot, templateCU, false);
            if (hasAdapterVariable) {
                insertInnerClass(iCompilationUnit, adapterInterfaceCU, false);
                insertInnerClass(iCU1, adapterImplCU1, true);
                insertInnerClass(iCU2, adapterImplCU2, true);
            }
            cleanImportDeclarations(iCompilationUnit);
            cleanImportDeclarations(iCU1);
            cleanImportDeclarations(iCU2);

        } else {
            // duplicate methods are in different packages
            log.info("duplicate methods are in different packages");
            return;
            /*
            saveCU(packageFragmentRoot, adapterInterfaceCU, false);
            saveCU(packageFragmentRoot, adapterImplCU1, false);
            saveCU(packageFragmentRoot, adapterImplCU2, false);
            saveCU(packageFragmentRoot, templateCU, false);
            */
        }

        saveMethod(iCU1, method1, cuImports1, true);
        saveMethod(iCU2, method2, cuImports2, true);

        packageFragmentRoot.getJavaProject().getProject().getWorkspace().save(true, null);
    }

    private ICompilationUnit saveCU(IPackageFragmentRoot packageFragmentRoot, CompilationUnit cu, boolean force)
            throws JavaModelException {
        String packageName = cu.getPackage().getName().getFullyQualifiedName();
        String name = ((TypeDeclaration) cu.types().get(0)).getName().getIdentifier() + ".java";
        String contents = cu.toString();
        IPackageFragment packageFragment = packageFragmentRoot.getPackageFragment(packageName);
        return packageFragment.createCompilationUnit(name, contents, force, null);
    }

    private void insertInnerClass(ICompilationUnit cu, CompilationUnit classCU, boolean inner) throws Exception {

        // creation of a Document
        Document document = new Document(cu.getSource());

        // creation of DOM/AST from a ICompilationUnit
        ASTParser parser = ASTParser.newParser(AST.JLS8);
        parser.setSource(cu);
        CompilationUnit astRoot = (CompilationUnit) parser.createAST(null);

        // creation of ASTRewrite
        ASTRewrite rewrite = ASTRewrite.create(astRoot.getAST());

        // description of the change
        ListRewrite classListRewrite;
        if (inner) {
            String typeDeclarationName = cu.getElementName().split("\\.")[0];
            TypeDeclaration typeDeclaration = null;
            List<TypeDeclaration> typeDeclarations = astRoot.types();
            for (TypeDeclaration curTypeDeclaration : typeDeclarations) {
                if (curTypeDeclaration.getName().getIdentifier().equals(typeDeclarationName)) {
                    typeDeclaration = curTypeDeclaration;
                }
            }

            if (typeDeclaration == null) {
                throw new IllegalStateException("cannot find target typeDeclaration: " + typeDeclarationName);
            }

            classListRewrite = rewrite.getListRewrite(typeDeclaration, TypeDeclaration.BODY_DECLARATIONS_PROPERTY);
        } else {
            classListRewrite = rewrite.getListRewrite(astRoot, CompilationUnit.TYPES_PROPERTY);
        }
        List<AbstractTypeDeclaration> types = classCU.types();
        for (AbstractTypeDeclaration type : types) {
            type.modifiers().clear();
            classListRewrite.insertLast(type, null);
        }
        ListRewrite imports = rewrite.getListRewrite(astRoot, IMPORTS_PROPERTY);
        List<ImportDeclaration> importDeclarations = classCU.imports();
        for (ImportDeclaration importDeclaration : importDeclarations) {
            imports.insertFirst(importDeclaration, null);
        }

        // computation of the text edits
        TextEdit edits = rewrite.rewriteAST(document, cu.getJavaProject().getOptions(true));

        // computation of the new source code
        edits.apply(document);
        String newSource = document.get();

        // update of the compilation unit
        cu.getBuffer().setContents(newSource);
        cu.getBuffer().save(null, true);
    }

    private void cleanImportDeclarations(ICompilationUnit cu) throws Exception {

        // creation of a Document
        Document document = new Document(cu.getSource());

        // creation of DOM/AST from a ICompilationUnit
        ASTParser parser = ASTParser.newParser(AST.JLS8);
        parser.setSource(cu);
        CompilationUnit astRoot = (CompilationUnit) parser.createAST(null);

        // creation of ASTRewrite
        ASTRewrite rewrite = ASTRewrite.create(astRoot.getAST());

        // description of the change
        ListRewrite imports = rewrite.getListRewrite(astRoot, IMPORTS_PROPERTY);
        Set<String> uniqueImportsSet = new HashSet<>();
        List<ImportDeclaration> importDeclarations = astRoot.imports();
        for (ImportDeclaration importDeclaration : importDeclarations) {
            String importDeclarationName = importDeclaration.getName().getFullyQualifiedName();
            if (importDeclaration.isOnDemand()) {
                importDeclarationName += ".*";
            }
            if (!uniqueImportsSet.add(importDeclarationName)) {
                imports.remove(importDeclaration, null);
            }
        }

        // computation of the text edits
        TextEdit edits = rewrite.rewriteAST(document, cu.getJavaProject().getOptions(true));

        // computation of the new source code
        edits.apply(document);
        String newSource = document.get();

        // update of the compilation unit
        cu.getBuffer().setContents(newSource);
        cu.getBuffer().save(null, true);

    }

    private void saveMethod(ICompilationUnit cu, MethodDeclaration method, List<Name> imports, boolean replace) throws Exception {

        // creation of a Document
        Document document = new Document(cu.getSource());

        // creation of DOM/AST from a ICompilationUnit
        ASTParser parser = ASTParser.newParser(AST.JLS8);
        parser.setSource(cu);
        CompilationUnit astRoot = (CompilationUnit) parser.createAST(null);

        IPath ignoreAnnotationPath = astRoot.getJavaElement().getPath();
        while(!ignoreAnnotationPath.lastSegment().equals("java")) {
        	ignoreAnnotationPath = ignoreAnnotationPath.removeLastSegments(1);
        }
        File ignoreAnnotationFile = new File("/Users/user/Desktop/JTestParametrizer-Benchmarks/gson" + ignoreAnnotationPath.addTrailingSeparator().append("ca/uwaterloo/jrefactoring/autogenerated/IgnoreClone.java").toString());
        if(!ignoreAnnotationFile.exists()) {
        	
        }

        // creation of ASTRewrite
        ASTRewrite rewrite = ASTRewrite.create(astRoot.getAST());

        Annotation ignoreCloneAnnotation = method.getAST().newMarkerAnnotation();
        Name ignoreCloneAnnotationName =  method.getAST().newSimpleName("ca");
        ignoreCloneAnnotationName = method.getAST().newQualifiedName(ignoreCloneAnnotationName, method.getAST().newSimpleName("uwaterloo"));
        ignoreCloneAnnotationName = method.getAST().newQualifiedName(ignoreCloneAnnotationName, method.getAST().newSimpleName("jrefactoring"));
        ignoreCloneAnnotationName = method.getAST().newQualifiedName(ignoreCloneAnnotationName, method.getAST().newSimpleName("autogenerated"));
        ignoreCloneAnnotationName = method.getAST().newQualifiedName(ignoreCloneAnnotationName, method.getAST().newSimpleName("IgnoreClone"));
        ignoreCloneAnnotation.setTypeName(ignoreCloneAnnotationName);

        // description of the change
        MethodVisitor methodVisitor = ASTNodeUtil.retrieveMethodDeclaration(astRoot, method);
        if (methodVisitor.getResult() != null) {
            if (replace) {
                ListRewrite methodDeclarations = rewrite.getListRewrite(methodVisitor.getTypeDeclaration(),
                        TypeDeclaration.BODY_DECLARATIONS_PROPERTY);
                methodDeclarations.replace(methodVisitor.getResult(), method, null);
                methodDeclarations.insertBefore(ignoreCloneAnnotation, method, null);
                ListRewrite cuImports = rewrite.getListRewrite(astRoot, IMPORTS_PROPERTY);
                for (Name importName : imports) {
                    ImportDeclaration importDeclaration = astRoot.getAST().newImportDeclaration();
                    importDeclaration.setName((Name) ASTNode.copySubtree(astRoot.getAST(), importName));
                    cuImports.insertFirst(importDeclaration, null);
                }
            } else {
                //log.info("repeated method found when the target method is not replaceable");
                throw new IllegalStateException("repeated method found when the target method is not replaceable");
            }
        } else {
            ListRewrite methodDeclarations = rewrite.getListRewrite(methodVisitor.getTypeDeclaration(),
                    TypeDeclaration.BODY_DECLARATIONS_PROPERTY);
            methodDeclarations.insertLast(ignoreCloneAnnotation, null);
            methodDeclarations.insertLast(method, null);
            ListRewrite cuImports = rewrite.getListRewrite(astRoot, IMPORTS_PROPERTY);
            for (Name importName : imports) {
                ImportDeclaration importDeclaration = astRoot.getAST().newImportDeclaration();
                importDeclaration.setName((Name) ASTNode.copySubtree(astRoot.getAST(), importName));
                cuImports.insertFirst(importDeclaration, null);
            }
        }

        // computation of the text edits
        TextEdit edits = rewrite.rewriteAST(document, cu.getJavaProject().getOptions(true));

        // computation of the new source code
        edits.apply(document);
        String newSource = document.get();

        // update of the compilation unit
        cu.getBuffer().setContents(newSource);
        cu.getBuffer().save(null, true);
    }

    @Override
    public String toString() {
        StringBuilder unrefactoredPairs = new StringBuilder();

        unrefactoredPairs.append("non-refactored node pairs: \n");
        for (NodePair nodePair : unrefactoredList) {
            unrefactoredPairs.append(nodePair.toString());
            unrefactoredPairs.append("\n");
        }

        return unrefactoredPairs.toString() + "\n" + (templateClass == null ? templateCUImports + "\n"
                + templateMethod.toString() : templateCU.toString()) + "\n"
                + adapterInterfaceCU.toString() + "\n"
                + adapterImplCU1.toString() + "\n"
                + adapterImplCU2.toString() + "\n"
                + cuImports1.toString() + "\n"
                + method1.toString() + "\n"
                + cuImports2.toString() + "\n"
                + method2.toString();
    }

}
