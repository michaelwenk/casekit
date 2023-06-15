package casekit.nmr.dbservice;

import casekit.nmr.model.DataSet;

import java.io.FileNotFoundException;
import java.util.List;

public class LOTUS {

    public static List<DataSet> getDataSetsWithShiftPredictionFromLOTUS(final String pathToLOTUS,
                                                                        final String[] nuclei) throws FileNotFoundException {
        final List<DataSet> dataSetList = COCONUT.getDataSetsWithShiftPredictionFromCOCONUT(pathToLOTUS, nuclei);
        for (final DataSet dataSet : dataSetList) {
            dataSet.getMeta()
                   .put("source", "lotus");
        }

        return dataSetList;
    }
}
