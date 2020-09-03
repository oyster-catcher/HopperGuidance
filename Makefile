# Standard Mac KSP install dir
KSP=/Users/${USER}/Library/Application\ Support/Steam/steamapps/common/Kerbal\ Space\ Program
VER=v0.2.4alpha

ASSEMBLYOPTS=-reference:${KSP}/KSP.app/Contents/Resources/Data/Managed/UnityEngine.CoreModule.dll -reference:${KSP}/KSP.app/Contents/Resources/Data/Managed/Assembly-CSharp.dll -reference:${KSP}/KSP.app/Contents/Resources/Data/Managed/UnityEngine.UI.dll -reference:${KSP}/KSP.app/Contents/Resources/Data/Managed/KSPAssets.dll -reference:${KSP}/KSP.app/Contents/Resources/Data/Managed/UnityEngine.PhysicsModule.dll -reference:${KSP}/KSP.app/Contents/Resources/Data/Managed/UnityEngine.InputLegacyModule.dll
SDKOPTS=-sdk:4.0

# 3rd Party - ALGLIB Library. See alglib.net
ALGLIB_SRCS=alglib/csharp/net-core/src/alglib_info.cs alglib/csharp/net-core/src/dataanalysis.cs alglib/csharp/net-core/src/interpolation.cs alglib/csharp/net-core/src/specialfunctions.cs alglib/csharp/net-core/src/alglibinternal.cs alglib/csharp/net-core/src/diffequations.cs alglib/csharp/net-core/src/linalg.cs alglib/csharp/net-core/src/statistics.cs alglib/csharp/net-core/src/alglibmisc.cs alglib/csharp/net-core/src/fasttransforms.cs alglib/csharp/net-core/src/optimization.cs alglib/csharp/net-core/src/ap.cs alglib/csharp/net-core/src/integration.cs alglib/csharp/net-core/src/Solvers.cs

.PHONY: all install

all: install Solve.exe HopperGuidance-${VER}.zip

install: GameData/HopperGuidance/Plugins/HopperGuidance.dll
	cp -r GameData ${KSP}
#	cp -r GameData ~/KSP_Cutdown

HopperGuidance-${VER}.zip: GameData/HopperGuidance/Plugins/HopperGuidance.dll
	rm -f HopperGuidance-${VER}.zip
	cd GameData; find HopperGuidance | zip -@ ../HopperGuidance-${VER}.zip

clean:
	rm -f HopperGuidance.dll *.exe *.zip

GameData/HopperGuidance/Plugins/HopperGuidance.dll: HopperGuidance.cs Solve.cs Trajectory.cs PID3d.cs GuiUtils.cs ConeUtils.cs
	mcs -define:UNITY ${ASSEMBLYOPTS} ${SDKOPTS} -target:library $^ ${ALGLIB_SRCS} -out:$@

Solve.exe : Solve.cs Trajectory.cs
	mcs ${ASSEMBLYOPTS} ${SDKOPTS} ${ALGLIB_SRCS} $^ -out:$@

ConeUtils.exe : ConeUtils.cs
	mcs ${ASSEMBLYOPTS} ${SDKOPTS} $^ -out:$@
