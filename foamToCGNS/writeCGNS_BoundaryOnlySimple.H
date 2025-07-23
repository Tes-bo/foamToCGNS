// 边界专用CGNS写入代码 - 简化版本（跳过场量处理）
{
#if defined _WRITE_CGNS_

    const pointField& points = mesh.points();
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    if ( access(cgnsDataPath.c_str(),W_OK) != 0 )
    {
        cerr << "请确保CGNS输出目录 " << cgnsDataPath << " 存在且可写。" << endl;
    }

    const std::string cgns_filename = cgnsDataPath + "/" + runTime.timeName()+"_boundaries_simple.cgns" ;

    Info << "CGNS输出文件(仅边界-简化版): " << cgns_filename << endl;

    // 创建CGNS文件
    CGNSOO::file cgnsFile( cgns_filename, CGNSOO::file::WRITE );
    const int cellDim = 2;  // 2D边界面
    const int physDim = 3;  // 3D物理空间
    CGNSOO::Base_t cgnsBase = cgnsFile.writeBase( "BoundaryBase", cellDim, physDim );

    // 使用SI单位
    cgnsBase.writeDataClass( CGNSOO::Dimensional );
    cgnsBase.writeSIUnits();

    // 稳态设置
    cgnsBase.writeSimulationType( CGNSOO::NonTimeAccurate );

    // ==============================================================
    // 步骤1: 收集所有边界点并创建映射
    // ==============================================================
    
    Info << "收集边界点..." << endl;
    
    labelHashSet boundaryPointsSet;
    
    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];
        if (pp.nPoints() > 0)
        {
            labelList patchPoints = pp.meshPoints();
            forAll(patchPoints, pointI)
            {
                boundaryPointsSet.insert(patchPoints[pointI]);
            }
        }
    }
    
    labelList boundaryPointsList = boundaryPointsSet.sortedToc();
    label nBoundaryPoints = boundaryPointsList.size();
    
    Info << "找到 " << nBoundaryPoints << " 个唯一的边界点" << endl;
    
    labelList pointMap(points.size(), -1);
    forAll(boundaryPointsList, i)
    {
        pointMap[boundaryPointsList[i]] = i;
    }
    
    // ==============================================================
    // 步骤2: 创建边界区域并写入坐标
    // ==============================================================
    
    int zoneNum = 0;
    std::ostringstream oss;
    oss << "BoundaryZone_" << zoneNum;
    std::string zoneName = oss.str();
    CGNSOO::ZoneType_t zoneType = CGNSOO::Unstructured;
    
    label nBoundaryFaces = 0;
    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];
        nBoundaryFaces += pp.size();
    }
    
    std::vector<int> vertexSize(1, nBoundaryPoints);
    std::vector<int> cellSize(1, nBoundaryFaces);
    std::vector<int> bndrySize;
    
    Info << "创建边界区域: " << vertexSize[0] << " 个节点, " << cellSize[0] << " 个面" << endl;

    CGNSOO::Zone_t cgnsZone = cgnsBase.writeZone( zoneName,
                              vertexSize, cellSize, bndrySize,
                              zoneType );

    // 写入边界坐标
    std::vector<double> xCoordinates(nBoundaryPoints);
    std::vector<double> yCoordinates(nBoundaryPoints);
    std::vector<double> zCoordinates(nBoundaryPoints);

    forAll(boundaryPointsList, i)
    {
        label originalPointI = boundaryPointsList[i];
        xCoordinates[i] = points[originalPointI].x();
        yCoordinates[i] = points[originalPointI].y();
        zCoordinates[i] = points[originalPointI].z();
    }

    CGNSOO::GridCoordinates_t cgnsGridCoord = cgnsZone.writeGridCoordinates();
    cgnsGridCoord.writeCoordinatesData( CGNSOO::GridCoordinates_t::CARTESIAN,
                        xCoordinates, yCoordinates, zCoordinates );

    Info << "坐标写入完成" << endl;

    // ==============================================================
    // 步骤3: 写入边界面连接性
    // ==============================================================
    
    Info << "处理边界面连接..." << endl;
    
    std::vector<int> quad_connectivity;
    std::vector<int> tri_connectivity;
    
    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];
        
        forAll(pp, faceI)
        {
            const face& f = pp[faceI];
            
            if (f.size() == 4)
            {
                quad_connectivity.push_back(pointMap[f[0]] + 1);
                quad_connectivity.push_back(pointMap[f[1]] + 1);
                quad_connectivity.push_back(pointMap[f[2]] + 1);
                quad_connectivity.push_back(pointMap[f[3]] + 1);
            }
            else if (f.size() == 3)
            {
                tri_connectivity.push_back(pointMap[f[0]] + 1);
                tri_connectivity.push_back(pointMap[f[1]] + 1);
                tri_connectivity.push_back(pointMap[f[2]] + 1);
            }
            else
            {
                for (label triI = 1; triI < f.size() - 1; triI++)
                {
                    tri_connectivity.push_back(pointMap[f[0]] + 1);
                    tri_connectivity.push_back(pointMap[f[triI]] + 1);
                    tri_connectivity.push_back(pointMap[f[triI + 1]] + 1);
                }
            }
        }
    }
    
    int startIndex = 1;
    
    if (!quad_connectivity.empty())
    {
        cgnsZone.writeElements(
            zoneName + std::string("_Quads"),
            CGNSOO::QUAD_4,
            startIndex,
            quad_connectivity.size()/4,
            0,
            quad_connectivity );
        startIndex += quad_connectivity.size()/4;
        Info << "写入 " << quad_connectivity.size()/4 << " 个四边形面" << endl;
    }
    
    if (!tri_connectivity.empty())
    {
        cgnsZone.writeElements(
            zoneName + std::string("_Tris"),
            CGNSOO::TRI_3,
            startIndex,
            tri_connectivity.size()/3,
            0,
            tri_connectivity );
        Info << "写入 " << tri_connectivity.size()/3 << " 个三角形面" << endl;
    }

    Info << "面连接性写入完成" << endl;

    // ==============================================================
    // 跳过场量处理 - 仅测试网格写入
    // ==============================================================

    Info << "跳过场量处理（简化版本）" << endl;

    // ==============================================================
    // 步骤4: 写入边界条件信息
    // ==============================================================

    Info << "处理来自 " << patches.size() << " 个补丁的边界条件" << endl;

    if (patches.size() > 0)
    {
        CGNSOO::ZoneBC_t cgnsZoneBC = cgnsZone.writeZoneBC();
        forAll(patches, patchI)
        {
            int noPatch = patchI+1;
            const polyPatch& pp = patches[patchI];

            if (pp.nPoints() == 0)
            {
                Info << noPatch << " --> 警告: 补丁名称: " << pp.name()
                    << " : 此补丁有零个点。跳过。" << endl;
            }
            else
            {
                labelList patch_index_points = pp.meshPoints();
                std::vector<int> bc_points_index;
                
                forAll(patch_index_points, pipI)
                {
                    label boundaryPointI = pointMap[patch_index_points[pipI]];
                    if (boundaryPointI >= 0)
                    {
                        bc_points_index.push_back(boundaryPointI + 1);
                    }
                }
                
                if (!bc_points_index.empty())
                {
                    word pt = pp.physicalType();
                    Info << noPatch << " : 写入边界补丁: " << pp.name() 
                         << " (" << bc_points_index.size() << " 个点)" << endl;
                    cgnsZoneBC.writeBC
                    (
                         pp.name(),
                         string_OpenFoam_to_CGNSBcType(pt),
                         CGNSOO::PointList,
                         bc_points_index
                    );
                }
            }
        }
    }

    Info << "边界专用CGNS导出完成（简化版本）。" << endl;

#endif
} 