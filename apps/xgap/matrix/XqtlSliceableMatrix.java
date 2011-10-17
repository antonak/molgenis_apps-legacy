package matrix;

import java.util.ArrayList;
import java.util.List;

import org.molgenis.framework.db.QueryRule;
import org.molgenis.framework.db.QueryRule.Operator;
import org.molgenis.matrix.MatrixException;
import org.molgenis.matrix.component.general.MatrixQueryRule;
import org.molgenis.matrix.component.interfaces.BasicMatrix;
import org.molgenis.matrix.component.interfaces.SliceableMatrix;


public class XqtlSliceableMatrix implements SliceableMatrix<String, String, Object>
{
	DataMatrixInstance wrappedMatrix;
	
	//keep track of names being used
	List<String> copiedRowNames;
	List<String> copiedColNames;
	
	//keep track of which indices are being used
	List<Integer> copiedRowIndices;
	List<Integer> copiedColIndices;
	
	// ????????
	int colOffset;
	int rowOffset;
	int colLimit;
	int rowLimit;
	
	public XqtlSliceableMatrix(DataMatrixInstance wrappedMatrix) throws MatrixException
	{
		this.wrappedMatrix = wrappedMatrix;
		this.reset();
	}

	@Override
	public List<String> getRowHeaders() throws MatrixException
	{
		return this.copiedRowNames;
	}

	@Override
	public List<String> getColHeaders() throws MatrixException
	{
		return this.copiedColNames;
	}

	@Override
	public List<Integer> getRowIndices() throws MatrixException
	{
		return this.copiedRowIndices;
	}

	@Override
	public List<Integer> getColIndices() throws MatrixException
	{
		return this.copiedColIndices;
	}
	
	public Object[] getRow(int index) throws Exception
	{
		if(!this.copiedRowIndices.contains(index))
		{
			throw new MatrixException("Index no longer in resultset");
		}
		
		Object[] wholeRow = wrappedMatrix.getRow(index);
		
		Object[] res = new Object[this.copiedRowIndices.size()];
		int resIndex = 0;
		for(Integer indexR : this.copiedRowIndices)
		{
			res[resIndex] = wholeRow[indexR];
			resIndex++;
		}

		return res;
	}
	
	public Object[] getCol(int index) throws Exception
	{
		if(!this.copiedColIndices.contains(index))
		{
			throw new MatrixException("Index no longer in resultset");
		}
		
		Object[] wholeCol = wrappedMatrix.getCol(index);
		
		Object[] res = new Object[this.copiedColIndices.size()];
		int resIndex = 0;
		for(Integer indexC : this.copiedColIndices)
		{
			res[resIndex] = wholeCol[indexC];
			resIndex++;
		}

		return res;
	}

	@Override
	public Object[][] getValues() throws MatrixException
	{
		try
		{
			return wrappedMatrix.getSubMatrix(this.copiedRowNames, this.copiedColNames).getElements();
		}
		catch(Exception e)
		{
			throw new MatrixException(e);
		}
	}

	@Override
	public Integer getColCount() throws MatrixException
	{
		return wrappedMatrix.getNumberOfCols();
	}

	@Override
	public Integer getRowCount() throws MatrixException
	{
		return wrappedMatrix.getNumberOfRows();
	}

	@Override
	public void refresh() throws MatrixException
	{
		this.reset();
	}

	@Override
	public List<? extends Object>[][] getValueLists() throws MatrixException
	{
		throw new UnsupportedOperationException();
	}

	@Override
	public SliceableMatrix<String, String, Object> slice(MatrixQueryRule rule) throws MatrixException
	{
		throw new UnsupportedOperationException();
	}

	@Override
	public SliceableMatrix<String, String, Object> sliceByColIndex(Operator operator, Integer index) throws Exception
	{
		throw new MatrixException("Unimplemented");
	}

	@Override
	public SliceableMatrix<String, String, Object> sliceByRowIndex(Operator operator, Integer index) throws Exception
	{
		throw new MatrixException("Unimplemented");
	}

	@Override
	public SliceableMatrix<String, String, Object> sliceByRowOffsetLimit(int limit, int offset) throws Exception
	{
		this.copiedRowNames = this.copiedRowNames.subList(offset, offset+limit);
		this.copiedRowIndices = this.copiedRowIndices.subList(offset, offset+limit);
		return this;
	}

	@Override
	public SliceableMatrix<String, String, Object> sliceByColOffsetLimit(int limit, int offset) throws Exception
	{
		throw new MatrixException("Unimplemented");
	}

	@Override
	public SliceableMatrix<String, String, Object> sliceByRowValues(int index, Operator operator, Object value)
			throws Exception
	{
		List<String> result = null;
		if (this.wrappedMatrix.getData().getValueType().equals("Decimal"))
		{
			double valueD = Double.parseDouble(value.toString());
			result = AbstractDataMatrixQueries.selectUsingDecimal(this.getRow(index), valueD, operator, this.copiedColNames);
		}
		else
		{
			String valueS = value.toString();
			result = AbstractDataMatrixQueries.selectUsingText(this.getRow(index), valueS, operator, this.copiedColNames);
		}
		
		if (result.size() == 0)
		{
			throw new Exception("No colnames in resultset, empty matrix!");
		}
		
		//update names
		this.copiedColNames = result;
		
		//update indices
		List<Integer> colIndices = new ArrayList<Integer>();
		for(String s : result)
		{
			colIndices.add(wrappedMatrix.getColIndexForName(s));
		}
		this.copiedColIndices = colIndices;
		
		return this;
	}

	@Override
	public SliceableMatrix<String, String, Object> sliceByRowValues(String row, Operator operator, Object value)
			throws Exception
	{
		throw new MatrixException("Unimplemented");
	}

	@Override
	public SliceableMatrix<String, String, Object> sliceByColValues(int index, Operator operator, Object value)
			throws Exception
	{
		throw new MatrixException("Unimplemented");
	}

	@Override
	public SliceableMatrix<String, String, Object> sliceByColValues(String col, Operator operator, Object value)
			throws Exception
	{
		throw new MatrixException("Unimplemented");
	}

	@Override
	public SliceableMatrix<String, String, Object> sliceByRowProperty(String property, Operator operator, Object value)
			throws MatrixException
	{
		throw new MatrixException("Unimplemented");
	}

	@Override
	public SliceableMatrix<String, String, Object> sliceByColProperty(String property, Operator operator, Object value)
			throws Exception
	{
		throw new MatrixException("Unimplemented");
	}

	@Override
	public SliceableMatrix<String, String, Object> sliceByColValueProperty(String col, String property,
			Operator operator, Object value) throws MatrixException
	{
		throw new MatrixException("Unimplemented");
	}

	@Override
	public SliceableMatrix<String, String, Object> sliceByColValueProperty(int colIndex, String property,
			Operator operator, Object value) throws MatrixException
	{
		throw new MatrixException("Unimplemented");
	}

	@Override
	public BasicMatrix<String, String, Object> getResult() throws Exception
	{
		throw new UnsupportedOperationException();
	}

	@Override
	public List<String> getRowPropertyNames()
	{
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public List<String> getColPropertyNames()
	{
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public List<String> getValuePropertyNames()
	{
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public void reset() throws MatrixException
	{
		this.copiedColNames = this.wrappedMatrix.getColNames();
		this.copiedRowNames = this.wrappedMatrix.getRowNames();
		
		for(int row = 0; row < copiedRowNames.size(); row++)
		{
			copiedRowIndices.add(row);
		}
		for(int col = 0; col < copiedColNames.size(); col++)
		{
			copiedColIndices.add(col);
		}
		
		this.colLimit = this.copiedColNames.size()-1;
		this.colOffset = 0;
		this.rowLimit = this.copiedRowNames.size()-1;
		this.rowOffset = 0;
		
	}

	@Override
	public int getRowLimit()
	{
		return this.rowLimit;
	}

	@Override
	public void setRowLimit(int rowLimit)
	{
		this.rowLimit = rowLimit;
	}

	@Override
	public int getRowOffset()
	{
		return this.rowOffset;
	}

	@Override
	public void setRowOffset(int rowOffset)
	{
		this.rowOffset = rowOffset;
	}

	@Override
	public int getColLimit()
	{
		return this.colLimit;
	}

	@Override
	public void setColLimit(int colLimit)
	{
		this.colLimit = colLimit;
	}

	@Override
	public int getColOffset()
	{
		return this.colOffset;
	}

	@Override
	public void setColOffset(int colOffset)
	{
		this.colOffset = colOffset;
	}

	@Override
	public SliceableMatrix<String, String, Object> sliceByRowValueProperty(String row, String property,
			Operator operator, Object value) throws MatrixException
	{
		throw new UnsupportedOperationException();
	}

	@Override
	public SliceableMatrix<String, String, Object> sliceByRowValueProperty(int rowIndex, String property,
			Operator operator, Object value) throws MatrixException
	{
		throw new UnsupportedOperationException();
	}


}
