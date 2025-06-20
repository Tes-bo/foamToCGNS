/*------------------------------------------------------------------------------
Copyright (C) 2004-2007 Hydro-Quebec

This file is part of CGNSOO

CGNSOO is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the
Free Software Foundation; either version 2, or (at your option) any
later version.

CGNSOO is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with CGNSOO  If not, write to the Free Software Foundation,
Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
------------------------------------------------------------------------------*/
#ifndef CGNSMLL_H
#define CGNSMLL_H

namespace CGNSOO
{

#include "cgnslib.h"

//Have to create the non-CG_ prefixed names, because the current library might
//be using e.g. CG_ZoneType_t
#if CG_BUILD_SCOPE
#define MassUnits_t CG_MassUnits_t
#define LengthUnits_t CG_LengthUnits_t
#define TimeUnits_t CG_TimeUnits_t
#define TemperatureUnits_t CG_TemperatureUnits_t
#define AngleUnits_t CG_AngleUnits_t
#define ElectricCurrentUnits_t CG_ElectricCurrentUnits_t
#define SubstanceAmountUnits_t CG_SubstanceAmountUnits_t
#define LuminousIntensityUnits_t CG_LuminousIntensityUnits_t
#define DataClass_t CG_DataClass_t
#define GridLocation_t CG_GridLocation_t
#define BCDataType_t CG_BCDataType_t
#define GridConnectivityType_t CG_GridConnectivityType_t
#define PointSetType_t CG_PointSetType_t
#define GoverningEquationsType_t CG_GoverningEquationsType_t
#define ModelType_t CG_ModelType_t
#define BCType_t CG_BCType_t
#define DataType_t CG_DataType_t
#define ElementType_t CG_ElementType_t
#define ZoneType_t CG_ZoneType_t
#define RigidGridMotionType_t CG_RigidGridMotionType_t
#define ArbitraryGridMotionType_t CG_ArbitraryGridMotionType_t
#define SimulationType_t CG_SimulationType_t
#define WallFunctionType_t CG_WallFunctionType_t
#define AreaType_t CG_AreaType_t
#define AverageInterfaceType_t CG_AverageInterfaceType_t

#define Unstructured CG_Unstructured
#define Structured CG_Structured

#define PointRange CG_PointRange
#define PointList CG_PointList
#define ElementRange CG_ElementRange
#define ElementList CG_ElementList
#define RealSingle CG_RealSingle
#define RealDouble CG_RealDouble
#define Integer CG_Integer

#define Kilogram CG_Kilogram
#define Meter CG_Meter
#define Second CG_Second
#define Kelvin CG_Kelvin
#define Degree CG_Degree
#define Character CG_Character

#define DataTypeUserDefined CG_DataTypeUserDefined
#define MIXED CG_MIXED

#define Dimensional CG_Dimensional
#define GoverningEquationsNull CG_GoverningEquationsNull
#define TwoEquation_JonesLaunder CG_TwoEquation_JonesLaunder

#define NSLaminar CG_NSLaminar
#define NSTurbulent CG_NSTurbulent
#define NSTurbulentIncompressible CG_NSTurbulentIncompressible

#define NonTimeAccurate CG_NonTimeAccurate
#define TimeAccurate CG_TimeAccurate

#define Abutting1to1 CG_Abutting1to1
#define PointListDonor CG_PointListDonor


#define BAR_2 CG_BAR_2
#define BAR_3 CG_BAR_3
#define BCTypeNull CG_BCTypeNull
#define ElementTypeNull CG_ElementTypeNull
#define ElementTypeUserDefined CG_ElementTypeUserDefined
#define HEXA_20 CG_HEXA_20
#define HEXA_27 CG_HEXA_27
#define HEXA_8 CG_HEXA_8
#define NGON_n CG_NGON_n
#define NODE CG_NODE
#define PENTA_15 CG_PENTA_15
#define PENTA_18 CG_PENTA_18
#define PENTA_6 CG_PENTA_6
#define PYRA_14 CG_PYRA_14
#define PYRA_5 CG_PYRA_5
#define QUAD_4 CG_QUAD_4
#define QUAD_8 CG_QUAD_8
#define QUAD_9 CG_QUAD_9
#define TETRA_10 CG_TETRA_10
#define TETRA_4 CG_TETRA_4
#define TRI_3 CG_TRI_3
#define TRI_6 CG_TRI_6

#define BCAxisymmetricWedge CG_BCAxisymmetricWedge
#define BCInflow CG_BCInflow
#define BCOutflow CG_BCOutflow
#define BCSymmetryPlane CG_BCSymmetryPlane
#define BCTypeUserDefined CG_BCTypeUserDefined
#define BCWall CG_BCWall
#define CellCenter CG_CellCenter
#define EdgeCenter CG_EdgeCenter
#define FaceCenter CG_FaceCenter
#define FamilySpecified CG_FamilySpecified
#define GridLocationNull CG_GridLocationNull
#define GridLocationUserDefined CG_GridLocationUserDefined
#define IFaceCenter CG_IFaceCenter
#define JFaceCenter CG_JFaceCenter
#define KFaceCenter CG_KFaceCenter
#define Vertex CG_Vertex
#endif

// Always define these:
#define MODE_READ CG_MODE_READ
#define MODE_WRITE CG_MODE_WRITE
#define MODE_MODIFY CG_MODE_MODIFY

#define UserDefined CG_UserDefined

// -------------------


typedef char cgnsstring[33]; //!< Representation of a text label in CGNS

inline void CGNSfree( char* p ) { cg_free(p); }
}

#endif
