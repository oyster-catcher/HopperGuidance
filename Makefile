# Standard Mac KSP install dir
KSP=/Users/${USER}/Library/Application\ Support/Steam/steamapps/common/Kerbal\ Space\ Program
VER=v0.2.4
MONODIR=/Library/Frameworks/Mono.framework/Versions/5.12.0

MCS=${MONODIR}/Commands/mcs


ASSEMBLYOPTS=-reference:${KSP}/KSP.app/Contents/Resources/Data/Managed/UnityEngine.CoreModule.dll -reference:${KSP}/KSP.app/Contents/Resources/Data/Managed/Assembly-CSharp.dll -reference:${KSP}/KSP.app/Contents/Resources/Data/Managed/UnityEngine.UI.dll -reference:${KSP}/KSP.app/Contents/Resources/Data/Managed/KSPAssets.dll -reference:${KSP}/KSP.app/Contents/Resources/Data/Managed/UnityEngine.PhysicsModule.dll -reference:${KSP}/KSP.app/Contents/Resources/Data/Managed/UnityEngine.InputLegacyModule.dll
SDKOPTS=-sdk:4.0

#ASSEMBLYOPTS=-reference:UnityEngine.CoreModule.dll -reference:UnityEngine.SharedInternalsModule.dll

# 3rd Party - ALGLIB Library. See alglib.net
ALGLIB_SRCS=alglib/csharp/net-core/src/alglib_info.cs alglib/csharp/net-core/src/dataanalysis.cs alglib/csharp/net-core/src/interpolation.cs alglib/csharp/net-core/src/specialfunctions.cs alglib/csharp/net-core/src/alglibinternal.cs alglib/csharp/net-core/src/diffequations.cs alglib/csharp/net-core/src/linalg.cs alglib/csharp/net-core/src/statistics.cs alglib/csharp/net-core/src/alglibmisc.cs alglib/csharp/net-core/src/fasttransforms.cs alglib/csharp/net-core/src/optimization.cs alglib/csharp/net-core/src/ap.cs alglib/csharp/net-core/src/integration.cs alglib/csharp/net-core/src/Solvers.cs

.PHONY: all install test

all: install SolveTest.exe VesselSim.exe HopperGuidance-${VER}.zip

install: GameData/HopperGuidance/Plugins/HopperGuidance.dll
	cp -r GameData ${KSP}
	cp -r GameData ~/KSP_Cutdown
	cp -r GameData ~/KSP_RO

HopperGuidance-${VER}.zip: GameData/HopperGuidance/Plugins/HopperGuidance.dll
	rm -f HopperGuidance-${VER}.zip
	cd GameData; find HopperGuidance | zip -@ ../HopperGuidance-${VER}.zip

clean:
	rm -f HopperGuidance.dll *.exe *.zip

GameData/HopperGuidance/Plugins/HopperGuidance.dll: HopperGuidance.cs Solve.cs Trajectory.cs PID3d.cs GuiUtils.cs ConeUtils.cs HGTypes.cs HGUtils.cs PDController.cs StockAeroUtil.cs
	mkdir -p GameData/HopperGuidance/Plugins
	${MCS} -define:UNITY ${ASSEMBLYOPTS} ${SDKOPTS} -target:library $^ ${ALGLIB_SRCS} -out:$@

SolveTest.exe : Solve.cs Trajectory.cs SolveTest.cs HGTypes.cs HGUtils.cs StockAeroUtil.cs
	${MCS} ${ASSEMBLYOPTS} ${SDKOPTS} ${ALGLIB_SRCS} $^ -out:$@

ConeUtils.exe : ConeUtils.cs
	${MCS} ${ASSEMBLYOPTS} ${SDKOPTS} $^ -out:$@

VesselSim.exe : VesselSim.cs HGTypes.cs HGUtils.cs Solve.cs Trajectory.cs PDController.cs PID3d.cs ConeUtils.cs StockAeroUtil.cs
	${MCS} ${ASSEMBLYOPTS} ${SDKOPTS} ${ALGLIB_SRCS} $^ -out:$@

Spin.exe : Spin.cs Quaternion.cs
	${MCS} ${ASSEMBLYOPTS} ${SDKOPTS} $^ -out:$@

RotBody.exe : RotBody.cs Quaternion.cs
	${MCS} ${ASSEMBLYOPTS} ${SDKOPTS} $^ -out:$@

test:
	./run_tests.sh
	./run_tests2.sh
	./run_tests3.sh
	./run_tests4.sh
