package org.molgenis.compute.workflowgenerator;

import org.molgenis.compute.ComputeApplication;
import org.molgenis.compute.ComputeFeature;
import org.molgenis.compute.ComputeProtocol;
import org.molgenis.compute.pipelinemodel.Pipeline;
import org.molgenis.compute.pipelinemodel.Script;
import org.molgenis.compute.pipelinemodel.Step;
import org.molgenis.compute.scriptserver.MCF;
import org.molgenis.compute.ui.ComputeAppPaths;
import org.molgenis.compute.ui.DatabaseUpdater;
import org.molgenis.compute.ui.DatabaseUpdaterGridGain;
import org.molgenis.compute.ui.DatabaseUpdaterSsh;
import org.molgenis.framework.db.Database;
import org.molgenis.framework.db.DatabaseException;
import org.molgenis.pheno.ObservedValue;
import org.molgenis.protocol.Workflow;
import org.molgenis.protocol.WorkflowElement;
import org.molgenis.protocol.WorkflowElementParameter;
import org.molgenis.util.HttpServletRequestTuple;
import org.molgenis.util.Tuple;

import javax.servlet.ServletContext;
import java.io.IOException;
import java.text.ParseException;
import java.text.SimpleDateFormat;
import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: georgebyelas
 * Date: 18/10/2011
 * Time: 10:43
 * To change this template use File | Settings | File Templates.
 */
public class WorkflowGenerator
{
    private boolean flagJustGenerate = false;

    private static final String LOG = "log";// reserved word for logging feature type used in ComputeFeature

    public static final String DATE_FORMAT_NOW = "yyyy-MM-dd-HH-mm-ss";
    private SimpleDateFormat sdf = null;

    //format to run pipeline in compute
    private Pipeline pipeline = null;
    private Step currentStep = null;
    private String strCurrentPipelineStep = "INITIAL";
    private int pipelineElementNumber = 0;

    private int stepNumber = 0;

    private List<ComputeFeature> allComputeFeatures = null;

    //compute
    private MCF mcf = null;
    private DatabaseUpdater updater = null;

    //map of all compute features/values
    private Hashtable<String, String> weavingValues = null;
    private HashMap<String, ComputeFeature> computeFeatures = new HashMap<String, ComputeFeature>();
    Hashtable<String, String> userValues = null;

    //whole workflow application
    private ComputeApplication wholeWorkflowApp = null;

    //some necessary values
    private Workflow target = null;
    private String applicationName = null;

    private ParameterWeaver weaver = new ParameterWeaver();

    private String remoteLocation = null;
    private boolean isToWriteLocally = false;
    private String localLocation = "/";

    public void processSingleWorksheet(Database db, Tuple request,
                                       Hashtable<String, String> userValues,
                                       Workflow workflow,
                                       String applicationName /* should be unique somehow */) throws Exception
    {
        this.userValues = userValues;
        this.target = workflow;
        this.applicationName = applicationName;

        if (!db.inTx())
            db.beginTx();

        if (mcf == null)
        {
            HttpServletRequestTuple req = (HttpServletRequestTuple) request;
            ServletContext servletContext = req.getRequest().getSession().getServletContext();
            mcf = (MCF) servletContext.getAttribute("MCF");

            createDatabaseUpdater(mcf);
        }

        System.out.println(">>> generate apps");

        //create new pipeline and set current step to null
        pipeline = new Pipeline();
        currentStep = null;
        stepNumber = 0;

        userValues = new Hashtable<String, String>();
        pipelineElementNumber = 0;

        //application for the whole workflow
        wholeWorkflowApp = new ComputeApplication();

        //get the chosen workflow
//        Workflow workflow = db.query(Workflow.class).find().get(0);
        wholeWorkflowApp.setProtocol(workflow);
        wholeWorkflowApp.setInterpreter("WorkflowInterpreter");

        //it would be nice to select compute features of only selected workflow
        allComputeFeatures = db.query(ComputeFeature.class).find();
        System.out.println("we have so many features: " + allComputeFeatures.size());

        System.out.println("workflow" + workflow.getName());

        //add few parameters
        wholeWorkflowApp.setTime(now());

        //set app name everywhere and add to database
        sdf = new SimpleDateFormat(DATE_FORMAT_NOW);
        wholeWorkflowApp.setName(applicationName);
        pipeline.setId(applicationName);
        weaver.setJobID(applicationName);
//        db.beginTx();
        db.add(wholeWorkflowApp);

        //process workflow elements
        List<WorkflowElement> workflowElements = db.query(WorkflowElement.class).equals(WorkflowElement.WORKFLOW, workflow.getId()).find();

        for (int i = 0; i < workflowElements.size(); i++)
        {
            WorkflowElement workflowElement = workflowElements.get(i);
            processWorkflowElement(db, request, workflowElement);
        }

        String logfile = weaver.getLogfilename();

        pipeline.setPipelinelogpath(logfile);

        db.commitTx();
        executePipeline(db, pipeline);
    }

    private void createDatabaseUpdater(MCF mcf)
    {
        if (mcf.getBasis().equalsIgnoreCase(MCF.GRID))
            updater = new DatabaseUpdaterGridGain(mcf);
        else if ((mcf.getBasis().equalsIgnoreCase(MCF.SSH)))
            updater = new DatabaseUpdaterSsh(mcf);
    }

    public Date now()
    {
        Calendar cal = Calendar.getInstance();
        return cal.getTime();
    }

    public void executePipeline(Database db, Pipeline pipeline)
    {
        if (mcf != null && !flagJustGenerate)
        {
            mcf.setPipeline(pipeline);

            if (!updater.isStarted())
            {
                updater.setSettings(10, 10);
                updater.setDatabase(db);
                updater.start();
            }
        }
        else
            System.out.println(pipeline.toString());

    }

    private void processWorkflowElement(Database db, Tuple request, WorkflowElement workflowElement)
            throws DatabaseException, ParseException, IOException
    {

        weavingValues = new Hashtable<String, String>();
        weavingValues.putAll(userValues);

        System.out.println(">>> workflow element: " + workflowElement.getName());

        //create complex features, which will be processed after simple features
        Vector<ComputeFeature> featuresToDerive = new Vector<ComputeFeature>();

        //features to iterate
        Vector<ComputeFeature> featuresToIterate = new Vector<ComputeFeature>();

        //get protocol and template
        ComputeProtocol protocol = db.findById(ComputeProtocol.class, workflowElement.getProtocol_Id());


        //process compute features
        for (ComputeFeature computeFeature : allComputeFeatures)
        {
            if (computeFeature.getIsUser())
                continue;
            else if (computeFeature.getIsDerived())
            {
                featuresToDerive.addElement(computeFeature);
            }
            else if (computeFeature.getIterateOver())
            {
                featuresToIterate.addElement(computeFeature);
            }
            else
            {
                weavingValues.put(computeFeature.getName(), computeFeature.getDefaultValue());
            }
        }


        //process workflow element parameters
        List<WorkflowElementParameter> workflowElementParameters = db.query(WorkflowElementParameter.class).
                equals(WorkflowElementParameter.WORKFLOWELEMENT, workflowElement.getId()).find();

        for (WorkflowElementParameter par : workflowElementParameters)
        {
            ComputeFeature feature = computeFeatures.get(par.getTarget_Name());
            weavingValues.put(par.getFeature_Name(), feature.getDefaultValue());
        }

        //still hardcoded fastQC iteration
        if (workflowElement.getName().equalsIgnoreCase("FastqcElement"))
        {
            for (int i = 1; i < 3; i++)
            {
                weavingValues.put(featuresToIterate.elementAt(0).getName(), i + "");
                generateComputeApplication(db, request, workflowElement, protocol, weavingValues, featuresToDerive);
            }
        }
        else
        {
            generateComputeApplication(db, request, workflowElement, protocol, weavingValues, featuresToDerive);
        }
    }

    private void generateComputeApplication(Database db, Tuple request,
                                            WorkflowElement workflowElement,
                                            ComputeProtocol protocol,
                                            Hashtable<String, String> weavingValues,
                                            Vector<ComputeFeature> featuresToDerive)
            throws IOException, DatabaseException, ParseException
    {
        ComputeApplication app = new ComputeApplication();
        app.setProtocol(protocol);
        app.setWorkflowElement(workflowElement);
        app.setTime(now());

        String appName = applicationName + "_" + workflowElement.getName() + "_" + pipelineElementNumber;
        app.setName(appName);
        System.out.println("---application---> " + appName);

        String protocolTemplate = protocol.getScriptTemplate();

//        System.out.println("--- template \n" + protocolTemplate);

        //weave complex features
        for (int i = 0; i < featuresToDerive.size(); i++)
        {
            ComputeFeature feature = featuresToDerive.elementAt(i);
            String featureName = feature.getName();
            String featureTemplate = feature.getDefaultValue();

            String featureValue = weaver.weaveFreemarker(featureTemplate, weavingValues);
//            System.out.println("complex-feature: " + featureName + " --> " + featureValue);
            weavingValues.put(featureName, featureValue);
        }

        String result = weaver.weaveFreemarker(protocolTemplate, weavingValues);
        app.setComputeScript(result);
        app.setInterpreter(protocol.getInterpreter());
        db.add(app);

        List<ComputeApplication> res = db.query(ComputeApplication.class).equals(ComputeApplication.NAME, app.getName()).find();
        if (res.size() != 1)
            throw new DatabaseException("ERROR while inserting into db");

        app = res.get(0);

        Set entries = weavingValues.entrySet();
        Iterator it = entries.iterator();

        //this is used for database update with ComputeAppPaths
        Vector<String> logpathfiles = new Vector<String>();

        while (it.hasNext())
        {
            Map.Entry entry = (Map.Entry) it.next();
            String name = (String) entry.getKey();
            String value = (String) entry.getValue();

            ObservedValue observedValue = new ObservedValue();
            observedValue.setValue(value);
            observedValue.setProtocolApplication(app);
            observedValue.setTarget(target.getId());
            ComputeFeature feature = computeFeatures.get(name);
            if (feature.getFeatureType().equalsIgnoreCase(LOG))
            {
                logpathfiles.addElement(value);
            }

            observedValue.setFeature(feature);
            db.add(observedValue);
        }

        pipelineElementNumber++;

        //create compute pipeline
        String scriptID = app.getName();
        weaver.setScriptID(scriptID);
        weaver.setActualCommand(result);

        weaver.setDefaults();

        if (protocol.getWalltime() != null)
            weaver.setWalltime(protocol.getWalltime());
        if (protocol.getClusterQueue() != null)
            weaver.setClusterQueue(protocol.getClusterQueue());
        if (protocol.getCores() != null)
            weaver.setCores(protocol.getCores() + "");
        if (protocol.getMemoryReq() != null)
            weaver.setMemoryReq(protocol.getMemoryReq() + "");

        //at some point of time can be added for the verification
        weaver.setVerificationCommand("\n");

        weaver.setDatasetLocation(remoteLocation);
        String scriptRemoteLocation = remoteLocation + "scripts/";

        String scriptFile = weaver.makeScript();

        String logfile = weaver.getLogfilename();
        pipeline.setPipelinelogpath(logfile);

        if(isToWriteLocally)
            weaver.writeToFile(localLocation + pipelineElementNumber + scriptID, scriptFile);

        //todo rewrite pipeline generation
        //look into proper choose of logfile and all path settings
        List<String> strPreviousWorkflowElements = workflowElement.getPreviousSteps_Name();

        Script pipelineScript = new Script(scriptID, scriptRemoteLocation, scriptFile.getBytes());

        if (protocol.getClusterQueue() != null)
            if (protocol.getClusterQueue().equalsIgnoreCase("short"))
                pipelineScript.setShort(true);

        if (strPreviousWorkflowElements.size() == 0)//script does not depend on other scripts
        {
            if (currentStep == null) //it is a first script in the pipeline
            {
                Step step = new Step(workflowElement.getName());
                step.setNumber(stepNumber);
                stepNumber++;
                currentStep = step;
                pipeline.addStep(step);
            }

            currentStep.addScript(pipelineScript);
        }
        else //scripts depends on previous scripts
        {
            String strPrevious = strPreviousWorkflowElements.get(0);

            if (!strPrevious.equalsIgnoreCase(strCurrentPipelineStep))
            {
                //Step step = new Step("step_" + app.getName());
                Step step = new Step(workflowElement.getName());
                step.setNumber(stepNumber);
                stepNumber++;
                currentStep = step;
                pipeline.addStep(step);
            }

            currentStep.addScript(pipelineScript);
            strCurrentPipelineStep = strPrevious;
        }

        //here ComputeAppPaths generation
        ComputeAppPaths appPaths = new ComputeAppPaths();
        appPaths.setApplication(app);
        appPaths.setErrpath(weaver.getErrfilename());
        appPaths.setOutpath(weaver.getOutfilename());
        appPaths.setExtralog(weaver.getExtralogfilename());

        if (logpathfiles.size() > 0)
            for (int iii = 0; iii < logpathfiles.size(); iii++)
                appPaths.addLogpath(logpathfiles.elementAt(iii));

        updater.addComputeAppPath(appPaths);
    }

    //root remote location should be set
    public void setRemoteLocation(String remoteLocation)
    {
        this.remoteLocation = remoteLocation;
    }

    public void setToWriteLocally(boolean toWriteLocally)
    {
        isToWriteLocally = toWriteLocally;
    }

    public void setLocalLocation(String localLocation)
    {
        this.localLocation = localLocation;
    }
}