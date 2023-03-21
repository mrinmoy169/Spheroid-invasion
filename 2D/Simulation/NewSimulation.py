
from cc3d import CompuCellSetup
        
from NewSimulationSteppables import NewSimulationSteppable
from NewSimulationSteppables import MitosisSteppable

CompuCellSetup.register_steppable(steppable=NewSimulationSteppable(frequency=1))
CompuCellSetup.register_steppable(steppable=MitosisSteppable(frequency=1))

CompuCellSetup.run()
