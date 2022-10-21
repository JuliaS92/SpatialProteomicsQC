import pytest
import pandas as pd
import domaps

#def test_weigh_yields_emptydf():
#    with pytest.raises(ValueError):
#        domaps.weigh_yields(pd.DataFrame())

@pytest.mark.parametrize("df, sets", [(pd.DataFrame([[1,0,1, 3,2,5, 3,15,20]], index=pd.Index(["a"],name="id"),
                                                     columns=pd.MultiIndex.from_arrays([["R","R","R","C","C","C","V","V","V"],
                                                                                        ["F1","F2","F3","F1","F2","F3","F1","F2","F3"]],
                                                                                       names=["Set", "Fraction"])),
                                        dict(ratio="Ra", counts="C", variability="V")),
                                       (pd.DataFrame([[1,0,1, 3,2,5, 3,15,20]], index=pd.Index(["a"],name="id"),
                                                     columns=pd.MultiIndex.from_arrays([["F1","F2","F3","F1","F2","F3","F1","F2","F3"]], names=["Fraction"])),
                                        dict(ratio="R", counts="C", variability="V")),
                                       (pd.DataFrame([[1,0,1, 3,2,5, 3,15,20]], index=pd.Index(["a"],name="id"),
                                                     columns=pd.MultiIndex.from_arrays([["R","R","R","C","C","C","V","V","V"],
                                                                                        ["F1","F2","F3","F1","F2","F3","F1","F2","F3"]],
                                                                                        names=["Set", "Fractions"])),
                                        dict(ratio="R", counts="C", variability="V"))])
def test_filter_SILAC_countvar_KeyErrors(df, sets):
    with pytest.raises(KeyError):
        domaps.filter_SILAC_countvar(df,sets)