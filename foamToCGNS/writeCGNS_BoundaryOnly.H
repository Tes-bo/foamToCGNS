// 边界专用CGNS写入代码 - 只导出边界网格和边界场量
// 边界专用CGNS写入代码 - 只导出边界网格和边界场量
{
    const polyBoundaryMesh& patches = mesh.boundaryMesh();
    const pointField& points = mesh.points();
    
    // 定义shouldExportPatch函数的逻辑，直接内联使用

    // 初始化输出目录
    std::string cgnsDataPath = runTime.path() + "/ConversionCGNS";
    if ( !isDir(cgnsDataPath) )
    {
        mkDir(cgnsDataPath);
    }

    if ( !isDir(cgnsDataPath) )
    {
        cerr << "请确保CGNS输出目录 " << cgnsDataPath << " 存在且可写。" << endl;
    }

    const std::string cgns_filename = cgnsDataPath + "/" + runTime.timeName()+"_boundaries.cgns" ;

    Info << "CGNS输出文件(仅边界): " << cgns_filename << endl;

    // 创建CGNS文件
    CGNSOO::file cgnsFile( cgns_filename, CGNSOO::file::WRITE );
    const int cellDim = 2;  // 2D边界面
    const int physDim = 3;  // 3D物理空间
    CGNSOO::Base_t cgnsBase = cgnsFile.writeBase( "BoundaryBase", cellDim, physDim );

#ifdef TURBOMACHINERY_WORKGROUP
    // 叶轮机械描述符
    cgnsBase.writeDescriptor( "Turbomachinery Workgroup - Boundary Only" );
#endif

    // 使用SI单位
    cgnsBase.writeDataClass( CGNSOO::Dimensional );
    cgnsBase.writeSIUnits();

    // 流动方程组
    CGNSOO::GoverningEquationsType_t equationSet = CGNSOO::GoverningEquationsNull;
#if defined CGNSTOFOAM_EXTRACT_FlowEquationSet
    if ( application == "simpleFoam" ) equationSet = CGNSOO::NSTurbulentIncompressible;
    else if ( application == "icoFoam" ) equationSet = CGNSOO::NSLaminar;
#endif

    CGNSOO::FlowEquationSet_t cgnsFlowEq = cgnsBase.writeFlowEquationSet(physDim);
    switch( equationSet )
    {
    case CGNSOO::NSTurbulentIncompressible:
        cgnsFlowEq.writeGoverningEquations( equationSet );
        cgnsFlowEq.writeTurbulenceModel( CGNSOO::TwoEquation_JonesLaunder );
        break;
    case CGNSOO::NSTurbulent:
        cgnsFlowEq.writeGoverningEquations( equationSet );
        cgnsFlowEq.writeTurbulenceModel( CGNSOO::TwoEquation_JonesLaunder );
        break;
    case CGNSOO::NSLaminar:
        cgnsFlowEq.writeGoverningEquations( equationSet );
        break;
    default:
        break;
    }

    // 稳态或非稳态
    cgnsBase.writeSimulationType( (steadyState) ? CGNSOO::NonTimeAccurate : CGNSOO::TimeAccurate );

    // ==============================================================
    // 步骤1: 收集所有边界点并创建映射 (跳过empty补丁)
    // ==============================================================
    
    Info << "收集边界点..." << endl;
    
    // 存储唯一边界点标签的集合
    labelHashSet boundaryPointsSet;
    
    // 收集指定的边界点
    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];
        word patchType = pp.type();
        
        // 检查是否应该导出此补丁
        bool shouldExport = true;
        
        // 跳过empty类型
        if (patchType == "empty")
        {
            shouldExport = false;
        }
        // 如果指定了补丁名称，只导出指定的补丁
        else if (selectedPatches.size() > 0)
        {
            bool found = false;
            forAll(selectedPatches, i)
            {
                if (pp.name() == selectedPatches[i])
                {
                    found = true;
                    break;
                }
            }
            shouldExport = found;
        }
        
        if (!shouldExport)
        {
            Info << "  跳过补丁: " << pp.name() << " (类型: " << patchType << ", 点数: " << pp.nPoints() << ", 面数: " << pp.size() << ")" << endl;
            continue;
        }
        
        Info << "  收集补丁: " << pp.name() << " (类型: " << patchType << ", 点数: " << pp.nPoints() << ", 面数: " << pp.size() << ")" << endl;
        
        if (pp.nPoints() > 0)
        {
            labelList patchPoints = pp.meshPoints();
            forAll(patchPoints, pointI)
            {
                boundaryPointsSet.insert(patchPoints[pointI]);
            }
        }
    }
    
    // 转换为排序列表以确保一致的顺序
    labelList boundaryPointsList = boundaryPointsSet.sortedToc();
    label nBoundaryPoints = boundaryPointsList.size();
    
    Info << "找到 " << nBoundaryPoints << " 个唯一的边界点" << endl;
    
    // 创建从原始点标签到新边界点索引的映射
    labelList pointMap(points.size(), -1);
    forAll(boundaryPointsList, i)
    {
        pointMap[boundaryPointsList[i]] = i;
    }
    
    // ==============================================================
    // 步骤2: 创建边界区域并写入坐标 (跳过empty补丁)
    // ==============================================================
    
    int zoneNum = 0;
    std::ostringstream oss;
    oss << "BoundaryZone_" << zoneNum;
    std::string zoneName = oss.str();
    CGNSOO::ZoneType_t zoneType = CGNSOO::Unstructured;
    
    // 计算指定边界面总数
    label nBoundaryFaces = 0;
    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];
        word patchType = pp.type();
        
        // 检查是否应该导出此补丁
        bool shouldExport = true;
        
        // 跳过empty类型
        if (patchType == "empty")
        {
            shouldExport = false;
        }
        // 如果指定了补丁名称，只导出指定的补丁
        else if (selectedPatches.size() > 0)
        {
            bool found = false;
            forAll(selectedPatches, i)
            {
                if (pp.name() == selectedPatches[i])
                {
                    found = true;
                    break;
                }
            }
            shouldExport = found;
        }
        
        if (!shouldExport)
        {
            continue;  // 跳过不需要导出的补丁
        }
        
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

    // ==============================================================
    // 步骤3: 写入边界面连接性 (跳过empty补丁)
    // ==============================================================
    
    Info << "处理边界面连接..." << endl;
    
    std::vector<int> quad_connectivity;
    std::vector<int> tri_connectivity;
    
    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];
        word patchType = pp.type();
        
        // 检查是否应该导出此补丁
        bool shouldExport = true;
        
        // 跳过empty类型
        if (patchType == "empty")
        {
            shouldExport = false;
        }
        // 如果指定了补丁名称，只导出指定的补丁
        else if (selectedPatches.size() > 0)
        {
            bool found = false;
            forAll(selectedPatches, i)
            {
                if (pp.name() == selectedPatches[i])
                {
                    found = true;
                    break;
                }
            }
            shouldExport = found;
        }
        
        if (!shouldExport)
        {
            continue;  // 跳过不需要导出的补丁
        }
        
        forAll(pp, faceI)
        {
            const face& f = pp[faceI];
            
            if (f.size() == 4)
            {
                // 四边形面
                quad_connectivity.push_back(pointMap[f[0]] + 1);
                quad_connectivity.push_back(pointMap[f[1]] + 1);
                quad_connectivity.push_back(pointMap[f[2]] + 1);
                quad_connectivity.push_back(pointMap[f[3]] + 1);
            }
            else if (f.size() == 3)
            {
                // 三角形面
                tri_connectivity.push_back(pointMap[f[0]] + 1);
                tri_connectivity.push_back(pointMap[f[1]] + 1);
                tri_connectivity.push_back(pointMap[f[2]] + 1);
            }
            else
            {
                // 将多边形转换为三角形
                for (label triI = 1; triI < f.size() - 1; triI++)
                {
                    tri_connectivity.push_back(pointMap[f[0]] + 1);
                    tri_connectivity.push_back(pointMap[f[triI]] + 1);
                    tri_connectivity.push_back(pointMap[f[triI + 1]] + 1);
                }
            }
        }
    }
    
    // 写入单元节
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

    // ==============================================================
    // 步骤4: 处理边界场数据 (跳过empty补丁)
    // ==============================================================

    Info << "处理边界场量..." << endl;

    HashSet<word> volScalarNames;
    HashSet<word> volVectorNames;
    HashSet<word> surfScalarNames;
    HashSet<word> surfVectorNames;
    HashSet<word> sprayScalarNames;
    HashSet<word> sprayVectorNames;

    #include "getFieldNames.H"

    PtrList<volScalarField>     volScalarFields;
    PtrList<volVectorField>     volVectorFields;

    // 读取体积标量场
    readFields<volScalarField, fvMesh> (mesh,
                        objects,
                        volScalarNames,
                        volScalarFields);

    // 读取体积矢量场
    readFields<volVectorField, fvMesh> (mesh,
                        objects,
                        volVectorNames,
                        volVectorFields);

    if ( volScalarFields.size() || volVectorFields.size() )
    {
        CGNSOO::FlowSolution_t cgnsSolution = cgnsZone.writeFlowSolution("BoundaryFlowSolution", CGNSOO::Vertex);

        // 处理标量场
        forAll(volScalarFields, fieldI)
        {
            const volScalarField& vsf = volScalarFields[fieldI];

            CGNSOO::Quantity_t e_cgns_qty = string_OpenFoam_to_CGNSField(vsf.name());
            if (e_cgns_qty == CGNSOO::NULL_DATA && !allow_userdefined_fields )
            {
                Info << "标量场 : " << vsf.name() << " 暂不支持" << endl;
                continue;
            }
            
            Info << "处理边界标量场 '" << vsf.name() << "'" << endl;

            // 为每个点提取边界值
            std::vector<double> boundaryValues(nBoundaryPoints, 0.0);
            std::vector<int> pointCount(nBoundaryPoints, 0);
            
            Info << "  初始化边界值数组，大小: " << nBoundaryPoints << endl;
            
            // 遍历补丁并提取场值
            forAll(patches, patchI)
            {
                const polyPatch& pp = patches[patchI];
                word physType = pp.physicalType();
                
                // 检查是否应该导出此补丁
                bool shouldExport = true;
                
                // 跳过empty类型
                if (physType == "empty")
                {
                    shouldExport = false;
                }
                // 如果指定了补丁名称，只导出指定的补丁
                else if (selectedPatches.size() > 0)
                {
                    bool found = false;
                    forAll(selectedPatches, i)
                    {
                        if (pp.name() == selectedPatches[i])
                        {
                            found = true;
                            break;
                        }
                    }
                    shouldExport = found;
                }
                
                if (!shouldExport)
                {
                    Info << "  跳过补丁: " << pp.name() << " (类型: " << pp.physicalType() << ")" << endl;
                    continue;  // 跳过不需要导出的补丁
                }
                
                Info << "  处理补丁 " << patchI << ": " << pp.name() << " (面数: " << pp.size() << ")" << endl;
                
                if (pp.size() == 0) continue;  // 跳过空补丁
                
                // 检查边界条件类型，跳过empty类型
                const fvPatchScalarField& patchField = vsf.boundaryField()[patchI];
                word bcType = patchField.type();
                
                if (bcType == "empty" || patchField.size() == 0)
                {
                    Info << "    跳过 " << bcType << " 类型的边界条件 (场值大小: " << patchField.size() << ")" << endl;
                    continue;
                }
                
                try 
                {
                    Info << "    获取补丁场成功，类型: " << bcType << ", 大小: " << patchField.size() << endl;
                    
                    forAll(pp, faceI)
                    {
                        if (faceI >= patchField.size())
                        {
                            Info << "    警告：面索引超出范围 " << faceI << " >= " << patchField.size() << endl;
                            continue;
                        }
                        
                        scalar faceValue = patchField[faceI];
                        const face& f = pp[faceI];
                        
                        forAll(f, pointI)
                        {
                            label globalPointI = f[pointI];
                            if (globalPointI >= pointMap.size())
                            {
                                Info << "    警告：点索引超出范围 " << globalPointI << " >= " << pointMap.size() << endl;
                                continue;
                            }
                            
                            label boundaryPointI = pointMap[globalPointI];
                            if (boundaryPointI >= 0 && boundaryPointI < nBoundaryPoints)
                            {
                                boundaryValues[boundaryPointI] += faceValue;
                                pointCount[boundaryPointI]++;
                            }
                        }
                    }
                    Info << "    补丁 " << pp.name() << " 处理完成" << endl;
                }
                catch (const std::exception& e)
                {
                    Info << "    错误：处理补丁 " << pp.name() << " 时发生异常: " << e.what() << endl;
                    continue;
                }
                catch (...)
                {
                    Info << "    错误：处理补丁 " << pp.name() << " 时发生未知异常" << endl;
                    continue;
                }
            }
            
            Info << "  所有补丁处理完成，开始平均值计算" << endl;
            
            // 平均点处的值
            for (label i = 0; i < nBoundaryPoints; i++)
            {
                if (pointCount[i] > 0)
                {
                    boundaryValues[i] /= pointCount[i];
                }
            }
            
            // 压力的特殊情况
            if (e_cgns_qty == CGNSOO::PRESSURE)
            {
                Info << "    按密度缩放压力 " << vsf.name() << " rho = " << rho_ << endl;
                for (label i = 0; i < nBoundaryPoints; i++)
                {
                    boundaryValues[i] *= rho_;
                }
            }

            std::string qty_s = (e_cgns_qty == CGNSOO::NULL_DATA)
                        ? std::string(vsf.name())
                        : (CGNSOO::QuantityEnumToString(e_cgns_qty));

            std::unique_ptr<std::vector<double> > p_data( new std::vector<double>(boundaryValues) );
            CGNSOO::DataArray_t field = cgnsSolution.writeField(qty_s, *p_data);

            if ( e_cgns_qty == CGNSOO::NULL_DATA )
            {
                dimensionSet dims = vsf.dimensions();
                std::vector<double> units(5);
                units[0] = dims[Foam::dimensionSet::MASS];
                units[1] = dims[Foam::dimensionSet::LENGTH];
                units[2] = dims[Foam::dimensionSet::TIME];
                units[3] = dims[Foam::dimensionSet::TEMPERATURE];
                units[4] = 0;
                field.writeDimensionalExponents(units);
            }
        }

        // 处理矢量场
        forAll(volVectorFields, fieldI)
        {
            const volVectorField& vvf = volVectorFields[fieldI];
            const std::string& offieldname = vvf.name();
            
            Info << "处理边界矢量场 '" << vvf.name() << "'" << endl;

            std::vector<double> boundaryValues_x(nBoundaryPoints, 0.0);
            std::vector<double> boundaryValues_y(nBoundaryPoints, 0.0);
            std::vector<double> boundaryValues_z(nBoundaryPoints, 0.0);
            std::vector<int> pointCount(nBoundaryPoints, 0);
            
            forAll(patches, patchI)
            {
                const polyPatch& pp = patches[patchI];
                word physType = pp.physicalType();
                
                // 检查是否应该导出此补丁
                bool shouldExport = true;
                
                // 跳过empty类型
                if (physType == "empty")
                {
                    shouldExport = false;
                }
                // 如果指定了补丁名称，只导出指定的补丁
                else if (selectedPatches.size() > 0)
                {
                    bool found = false;
                    forAll(selectedPatches, i)
                    {
                        if (pp.name() == selectedPatches[i])
                        {
                            found = true;
                            break;
                        }
                    }
                    shouldExport = found;
                }
                
                if (!shouldExport)
                {
                    continue;  // 跳过不需要导出的补丁
                }
                
                const fvPatchVectorField& patchField = vvf.boundaryField()[patchI];
                
                // 检查边界条件类型，跳过empty类型
                word bcType = patchField.type();
                
                if (bcType == "empty" || patchField.size() == 0)
                {
                    Info << "    跳过矢量场 " << bcType << " 类型的边界条件 (场值大小: " << patchField.size() << ")" << endl;
                    continue;
                }
                
                forAll(pp, faceI)
                {
                    if (faceI >= patchField.size())
                    {
                        Info << "    警告：矢量场面索引超出范围 " << faceI << " >= " << patchField.size() << endl;
                        continue;
                    }
                    
                    Foam::vector faceValue = patchField[faceI];
                    const face& f = pp[faceI];
                    
                    forAll(f, pointI)
                    {
                        label globalPointI = f[pointI];
                        label boundaryPointI = pointMap[globalPointI];
                        if (boundaryPointI >= 0)
                        {
                            boundaryValues_x[boundaryPointI] += faceValue.x();
                            boundaryValues_y[boundaryPointI] += faceValue.y();
                            boundaryValues_z[boundaryPointI] += faceValue.z();
                            pointCount[boundaryPointI]++;
                        }
                    }
                }
            }
            
            for (label i = 0; i < nBoundaryPoints; i++)
            {
                if (pointCount[i] > 0)
                {
                    boundaryValues_x[i] /= pointCount[i];
                    boundaryValues_y[i] /= pointCount[i];
                    boundaryValues_z[i] /= pointCount[i];
                }
            }

            CGNSOO::Quantity_t e_cgns_qty = string_OpenFoam_to_CGNSField(offieldname+"_X");
            if ( e_cgns_qty == CGNSOO::NULL_DATA && !allow_userdefined_fields )
            {
                Info << "矢量场 '" << vvf.name() << "' 暂不支持" << endl;
                continue;
            }

            CGNSOO::Quantity_t e_cgns_qty_x = string_OpenFoam_to_CGNSField(offieldname+"_X");
            CGNSOO::Quantity_t e_cgns_qty_y = string_OpenFoam_to_CGNSField(offieldname+"_Y");
            CGNSOO::Quantity_t e_cgns_qty_z = string_OpenFoam_to_CGNSField(offieldname+"_Z");
            std::string qty_x_s = (e_cgns_qty_x == CGNSOO::NULL_DATA)
                        ? vvf.name() + "_X"
                        : CGNSOO::QuantityEnumToString(e_cgns_qty_x);
            std::string qty_y_s = (e_cgns_qty_y == CGNSOO::NULL_DATA)
                        ? vvf.name() + "_Y"
                        : CGNSOO::QuantityEnumToString(e_cgns_qty_y);
            std::string qty_z_s = (e_cgns_qty_z == CGNSOO::NULL_DATA)
                        ? vvf.name() + "_Z"
                        : CGNSOO::QuantityEnumToString(e_cgns_qty_z);

            std::unique_ptr<std::vector<double> > p_data_x( new std::vector<double>(boundaryValues_x) );
            std::unique_ptr<std::vector<double> > p_data_y( new std::vector<double>(boundaryValues_y) );
            std::unique_ptr<std::vector<double> > p_data_z( new std::vector<double>(boundaryValues_z) );

            CGNSOO::DataArray_t fx = cgnsSolution.writeField(qty_x_s, *p_data_x);
            CGNSOO::DataArray_t fy = cgnsSolution.writeField(qty_y_s, *p_data_y);
            CGNSOO::DataArray_t fz = cgnsSolution.writeField(qty_z_s, *p_data_z);

            if ( e_cgns_qty == CGNSOO::NULL_DATA )
            {
                dimensionSet dims = vvf.dimensions();
                std::vector<double> units(5);
                units[0] = dims[Foam::dimensionSet::MASS];
                units[1] = dims[Foam::dimensionSet::LENGTH];
                units[2] = dims[Foam::dimensionSet::TIME];
                units[3] = dims[Foam::dimensionSet::TEMPERATURE];
                units[4] = 0;
                fx.writeDimensionalExponents( units );
                fy.writeDimensionalExponents( units );
                fz.writeDimensionalExponents( units );
            }
        }
    }

    // ==============================================================
    // 步骤5: 写入边界条件信息 (跳过empty补丁)
    // ==============================================================

    Info << "处理来自 " << patches.size() << " 个补丁的边界条件" << endl;

    if (patches.size() > 0)
    {
        CGNSOO::ZoneBC_t cgnsZoneBC = cgnsZone.writeZoneBC();
        forAll(patches, patchI)
        {
            int noPatch = patchI+1;
            const polyPatch& pp = patches[patchI];
            word physType = pp.physicalType();
            
            // 检查是否应该导出此补丁
            bool shouldExport = true;
            
            // 跳过empty类型
            if (physType == "empty")
            {
                shouldExport = false;
            }
            // 如果指定了补丁名称，只导出指定的补丁
            else if (selectedPatches.size() > 0)
            {
                bool found = false;
                forAll(selectedPatches, i)
                {
                    if (pp.name() == selectedPatches[i])
                    {
                        found = true;
                        break;
                    }
                }
                shouldExport = found;
            }
            
            if (!shouldExport)
            {
                Info << noPatch << " --> 跳过补丁: " << pp.name() << " (类型: " << pp.physicalType() << ")" << endl;
                continue;  // 跳过不需要导出的补丁
            }

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

    Info << "边界专用CGNS导出完成。" << endl;

}