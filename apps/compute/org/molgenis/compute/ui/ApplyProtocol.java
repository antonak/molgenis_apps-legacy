
package org.molgenis.compute.ui;

import java.text.ParseException;
import java.util.ArrayList;
import java.util.List;

import org.molgenis.framework.db.Database;
import org.molgenis.framework.db.DatabaseException;
import org.molgenis.framework.ui.EasyPluginController;
import org.molgenis.framework.ui.ScreenController;
import org.molgenis.ngs.NgsSample;
import org.molgenis.pheno.Measurement;
import org.molgenis.pheno.ObservedValue;
import org.molgenis.protocol.Protocol;
import org.molgenis.protocol.ProtocolApplication;
import org.molgenis.util.Tuple;

/**
 * ApplyProtocolController takes care of all user requests and application logic.
 *
 * <li>Each user request is handled by its own method based action=methodName. 
 * <li> MOLGENIS takes care of db.commits and catches exceptions to show to the user
 * <li>ApplyProtocolModel holds application state and business logic on top of domain model. Get it via this.getModel()/setModel(..)
 * <li>ApplyProtocolView holds the template to show the layout. Get/set it via this.getView()/setView(..).
 */
public class ApplyProtocol extends EasyPluginController<ApplyProtocolModel>
{	
	boolean add = true;
	//this is a model + controller ;-)
	Protocol protocol;
	List<Measurement> measurements;
	List<NgsSample> samples;
	ProtocolApplication application;
	List<ObservedValue> values;
	
	public ApplyProtocol(String name, ScreenController<?> parent)
	{
		super(name, null, parent);
		this.setModel(new ApplyProtocolModel(this));
		//set default view to StartView
		this.setView(new ApplyProtocolStartView(this)); 
	}
	
	@Override
	public void reload(Database db) throws Exception
	{	

	}
	
	public void createProtocolApplication(Database db, Tuple request) throws DatabaseException, ParseException
	{
		//find protocol
		protocol = Protocol.findById(db, request.getInt("Protocol"));
		
		//find the measurement attached to protocol
		measurements = db.query(Measurement.class).in(Measurement.ID, protocol.getFeatures_Id()).find();
		
		//find the samples
		samples = db.query(NgsSample.class).in(NgsSample.ID, request.getList("Samples")).find();
		
		//create a new protocol application using the selected protocol and samples.
		application = new ProtocolApplication();
		application.setName(protocol.getName() + System.currentTimeMillis());
		application.setProtocol(protocol);
		//db.add(application);
		
		//create empty observedValue for each of the features (but don't save them unless edited.
		values = new ArrayList<ObservedValue>();
		for(Measurement m: measurements)
		{
			for(NgsSample s: samples)
			{
				ObservedValue v = new ObservedValue();
				v.setTarget(s);
				v.setFeature(m);
				v.setProtocolApplication(application);
				
				values.add(v);
			}
		}
		
		this.setView(new ApplyProtocolEditView(this));
	}
	
	public void restart(Database db, Tuple request)
	{
		protocol = null;
		measurements = null;
		samples = null;
		application = null;
		values = null;
		this.setView(new ApplyProtocolStartView(this)); 
	}
	
	public void cancel(Database db, Tuple request)
	{
		this.restart(db,request);
	}
	
	public void saveProtocolApplication(Database db, Tuple request) throws Exception
	{
		//update application
		this.application.set(request);
		
		//update values
		//value_<measurement>_<sample>
		for(ObservedValue v: values)
		{
			v.setValue(request.getString("value_"+v.getTarget_Id()+"_"+v.getFeature_Id()));
		}
		
		if(add)
		{
			db.add(application);
			db.add(values);
			add = false;
		}
		else
		{
			db.update(application);
			db.update(application);
		}
		
		this.setView(new ApplyProtocolStartView(this));
	}

}