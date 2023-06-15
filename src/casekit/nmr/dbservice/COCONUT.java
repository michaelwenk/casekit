package casekit.nmr.dbservice;

import casekit.nmr.model.*;
import casekit.nmr.utils.Utils;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.io.iterator.IteratingSDFReader;
import org.openscience.cdk.silent.SilentChemObjectBuilder;

import java.io.FileNotFoundException;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.List;

public class COCONUT {

    public static List<DataSet> getDataSetsWithShiftPredictionFromCOCONUT(final String pathToCOCONUT,
                                                                          final String[] nuclei) throws FileNotFoundException {
        final List<DataSet> dataSetList = new ArrayList<>();
        final IteratingSDFReader iterator = new IteratingSDFReader(new FileReader(pathToCOCONUT),
                                                                   SilentChemObjectBuilder.getInstance());
        IAtomContainer structure;
        DataSet dataSet;
        Spectrum spectrum;
        Assignment assignment;
        String[] split, split2;
        String spectrumPropertyString, multiplicity;
        double calcShift;
        List<Integer> closestSignalList;
        int atomIndex;

        while (iterator.hasNext()) {
            structure = iterator.next();
            try {
                dataSet = Utils.atomContainerToDataSet(structure);

                for (final String nucleus : nuclei) {
                    final String atomType = casekit.nmr.utils.Utils.getAtomTypeFromNucleus(nucleus);
                    final List<Integer> atomIndices = Utils.getAtomTypeIndicesByElement(structure, atomType);
                    spectrumPropertyString = structure.getProperty("Predicted "
                                                                           + nucleus
                                                                           + " shifts", String.class);
                    if (spectrumPropertyString
                            == null) {
                        spectrumPropertyString = structure.getProperty("Predicted_"
                                                                               + nucleus
                                                                               + "_shifts", String.class);
                    }
                    if (spectrumPropertyString
                            == null) {
                        continue;
                    }

                    spectrumPropertyString = spectrumPropertyString.replaceAll("[\\n\\r]", ";");
                    split = spectrumPropertyString.split(";");
                    spectrum = new Spectrum();
                    spectrum.setNuclei(new String[]{nucleus});
                    spectrum.setSignals(new ArrayList<>());
                    assignment = new Assignment();
                    assignment.setNuclei(spectrum.getNuclei());
                    assignment.initAssignments(spectrum.getSignalCount());
                    for (int i = 0; i
                            < split.length; i++) {
                        atomIndex = atomIndices.get(i);
                        multiplicity = Utils.getMultiplicityFromProtonsCount(structure.getAtom(atomIndex)
                                                                                      .getImplicitHydrogenCount());
                        if (multiplicity
                                == null) {
                            continue;
                        }
                        multiplicity = multiplicity.toLowerCase();
                        split2 = split[i].split("\\s+");
                        calcShift = Double.parseDouble(split2[1]);
                        // add assignment (at first here because of search for already existing equivalent signals)
                        // just to be sure that we take the right signal if equivalences are present
                        closestSignalList = spectrum.pickByClosestShift(calcShift, 0, 0.0);
                        closestSignalList.retainAll(spectrum.pickByMultiplicity(multiplicity));
                        if (closestSignalList.isEmpty()) {
                            assignment.addAssignment(0, new int[]{atomIndex});
                        } else {
                            assignment.addAssignmentEquivalence(0, closestSignalList.get(0), atomIndex);
                        }
                        // add signal
                        spectrum.addSignal(
                                new Signal(new String[]{nucleus}, new Double[]{calcShift}, multiplicity, "signal", null,
                                           1, 0, null, null));
                    }

                    // if no spectrum could be built or the number of signals in spectrum is different than the atom number in molecule
                    if (Utils.getDifferenceSpectrumSizeAndMolecularFormulaCount(spectrum,
                                                                                Utils.getMolecularFormulaFromString(
                                                                                        dataSet.getMeta()
                                                                                               .get("mf")), 0)
                            != 0) {
                        continue;
                    }
                    dataSet.setSpectrum(new SpectrumCompact(spectrum));
                    dataSet.setAssignment(assignment);

                    dataSetList.add(dataSet.buildClone());
                }
            } catch (final CDKException e) {
                System.out.println("--> could not create a proper dataset for ID: "
                                           + structure.getProperty("ID"));
                e.printStackTrace();
            }
        }

        return dataSetList;
    }

}
