package casekit.nmr.elucidation.model;

import lombok.AllArgsConstructor;
import lombok.Getter;
import lombok.NoArgsConstructor;
import lombok.Setter;

@NoArgsConstructor
@AllArgsConstructor
@Getter
@Setter
public class ElucidationOptions {

    // PyLSD options
    private String[] filterPaths;
    private String pathToNeighborsFiles;
    private String pathToFragmentFiles;
    private boolean allowHeteroHeteroBonds;
    private boolean useElim;
    private int elimP1;
    private int elimP2;
}
