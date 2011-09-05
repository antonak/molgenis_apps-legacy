package regressiontest.archives;

import java.io.File;

import org.molgenis.framework.db.Database;
import org.molgenis.framework.db.DatabaseException;
import org.molgenis.organization.Investigation;
import org.molgenis.util.TarGz;
import org.molgenis.xgap.xqtlworkbench.ResetXgapDb;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import filehandling.storage.StorageHandler;

import plugins.archiveexportimport.XgapExcelExport;

import app.servlet.MolgenisServlet;

/**
 * Test data matrix import and export across all backends, all retrieval functions,
 * data types, and most dimensions, transpositions, sparsities, and text length variation.
 * 
 * To be used in xQTL automated test cases
 *
 */
public class TestArchives {

	@BeforeClass
	public void setup() throws Exception {
		
		Database db = new MolgenisServlet().getDatabase();
		StorageHandler sh = new StorageHandler(db);
		
		//assert db is empty
		Assert.assertFalse(sh.hasFileStorage(false));
		try{
			db.find(Investigation.class).get(0);
			Assert.fail("DatabaseException expected");
		}catch(DatabaseException expected){
			//DatabaseException was thrown
		}
		
		//setup database
		String report = ResetXgapDb.reset(db, true);
		Assert.assertTrue(report.endsWith("SUCCESS"));
		
		//setup file storage
		String path = "./tmp_archives_test_data";
		sh.setFileStorage(path);
		sh.validateFileStorage();
		Assert.assertTrue(sh.hasValidFileStorage());
	}

	@Test
	public void exportExcelArchive() throws Exception{
		System.out.println("** exportExcelArchive **");
		Thread.sleep(1000);
		Database db = new MolgenisServlet().getDatabase();
		File tmpDir = new File(System.getProperty("java.io.tmpdir") + File.separator
				+ "everyinvestigation" + "_export_" + System.nanoTime());
		tmpDir.mkdir();
		new XgapExcelExport(tmpDir, db);
		File tarFile = TarGz.tarDir(tmpDir);
		System.out.println("** Tarred into: "+tarFile.getAbsolutePath()+" **");
	}
	
	@Test
	public void importExcelArchive() throws Exception{
		Database db = new MolgenisServlet().getDatabase();
		
	}
	
}