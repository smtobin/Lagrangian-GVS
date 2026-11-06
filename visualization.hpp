#pragma once

#include "common.hpp"
#include "robot.hpp"

#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkPolyDataNormals.h>
#include <vtkTriangle.h>
#include <vtkActor.h>
#include <vtkCamera.h>
#include <vtkCubeSource.h>
#include <vtkSphereSource.h>
#include <vtkPlaneSource.h>
#include <vtkHDRReader.h>
#include <vtkImageFlip.h>
#include <vtkImageReader2.h>
#include <vtkImageReader2Factory.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkLight.h>
#include <vtkNamedColors.h>
#include <vtkOpenGLRenderer.h>
#include <vtkOpenGLTexture.h>
#include <vtkPBRIrradianceTexture.h>
#include <vtkPNGReader.h>
#include <vtkPNGWriter.h>
#include <vtkWindowToImageFilter.h>
#include <vtkPolyDataMapper.h>
#include <vtkPolyDataTangents.h>
#include <vtkProperty.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkSkybox.h>
#include <vtkTexture.h>
#include <vtkTriangleFilter.h>

#include <vtkSequencePass.h>
#include <vtkShadowMapBakerPass.h>
#include <vtkShadowMapPass.h>
#include <vtkCameraPass.h>
#include <vtkRenderPassCollection.h>
#include <vtkRenderStepsPass.h>
#include <vtkToneMappingPass.h>
#include <vtkLightsPass.h>
#include <vtkOpaquePass.h>
#include <vtkTranslucentPass.h>
#include <vtkImageShiftScale.h>
#include <vtkDepthPeelingPass.h>

#include <vtkCoordinate.h>

/** Creates render geometry for continuum robot */
struct ContinuumRobotGraphicsObject
{
    /** Number of sample points along continuum robot */
    unsigned num_samples;

    /** Pointer to continuum robot being visualized */
    ContinuumRobot* robot;

private:
    /** Cache vector of cross section points */
    std::vector<Vec3r> _cross_section_points;

    /** Integration points along robot - evenly spaced from 0 to L */
    std::vector<Real> _s_ips;

    /** Whether or not to draw end caps */
    bool _draw_end_caps;

    /** Resolution of generated geometry */
    int _cap_resolution;
    int _tubular_resolution;

    /** VTK data */
    vtkSmartPointer<vtkActor> _vtk_actor;
    vtkSmartPointer<vtkPolyData> _vtk_poly_data;
    vtkSmartPointer<vtkPolyDataNormals> _vtk_poly_data_normals;
    vtkSmartPointer<vtkPolyDataMapper> _vtk_poly_data_mapper;

public:
    ContinuumRobotGraphicsObject(ContinuumRobot* robot_, unsigned num_samples_)
        : robot(robot_)
        , num_samples(num_samples_)
        , _s_ips(num_samples)
        , _draw_end_caps(false)
        , _cap_resolution(10)
        , _tubular_resolution(15)
    {
        for (unsigned i = 0; i < num_samples; i++)
        {
            _s_ips[i] = Real(i) / (num_samples - 1) * robot->_L;
        }

        _vtk_poly_data = vtkSmartPointer<vtkPolyData>::New();
        
        vtkNew<vtkPolyDataMapper> data_mapper;
        // smooth normals
        vtkNew<vtkPolyDataNormals> normal_generator;
        normal_generator->SetInputData(_vtk_poly_data);
        normal_generator->SetFeatureAngle(30.0);
        normal_generator->SplittingOff();
        normal_generator->ConsistencyOn();
        normal_generator->ComputePointNormalsOn();
        normal_generator->ComputeCellNormalsOff();
        normal_generator->Update();

        data_mapper->SetInputConnection(normal_generator->GetOutputPort());

        _generateInitialPolyData();

        _vtk_actor = vtkSmartPointer<vtkActor>::New();
        _vtk_actor->SetMapper(data_mapper);
    }

    vtkSmartPointer<vtkActor> actor() { return _vtk_actor; }

    void _generateInitialPolyData()
    {
        // integrate initial robot state
        std::vector<Mat4r> Ts = robot->integratePoseAlongLength(_s_ips, robot->_q);

        vtkNew<vtkPoints> points;

        int _cap_resolution = 10;
        int _tubular_resolution = 10;

        // base cap
        Mat3r R_beg = Ts.front().block<3,3>(0,0);
        Vec3r p_beg = Ts.front().block<3,1>(0,3);
        if (_draw_end_caps)
        {
            for (int i = _cap_resolution-1; i >= 0; i--)
            {
                Real phi = (M_PI / 2.0) * i / _cap_resolution;

                Real ring_radius = robot->_r * std::cos(phi);
                Real z = robot->_r * std::sin(phi);

                for (int j = 0; j < _tubular_resolution; j++)
                {
                    Real theta = 2.0 * M_PI * j / _tubular_resolution;

                    Vec3r p(
                        ring_radius * std::cos(theta),
                        ring_radius * std::sin(theta),
                        -z);

                    Vec3r transformed = R_beg * p + p_beg;
                    points->InsertNextPoint(transformed.data());
                }
            }
        }

        // create cross-section points for a circular cross-section
        _cross_section_points.resize(_tubular_resolution);
        for (int i = 0; i < _tubular_resolution; i++)
        {
            const Real angle = i * 2 * M_PI / _tubular_resolution;
            _cross_section_points[i] = Vec3r( robot->_r * std::cos(angle), robot->_r * std::sin(angle), 0);
        }

        // insert cross-section vertices
        for (int si = 0; si < num_samples; si++)
        {
            Vec3r position = Ts[si].block<3,1>(0,3);
            Mat3r orientation = Ts[si].block<3,3>(0,0);
            // for (unsigned pi = 0; pi < cross_section_points.size(); pi++)
            for (const auto& p : _cross_section_points)
            {
                // transform each point in local cross section to global frame
                const Vec3r transformed = orientation * p + position;
                points->InsertNextPoint( transformed[0], transformed[1], transformed[2] );
            }
        }

        // end cap
        Mat3r R_end = Ts.back().block<3,3>(0,0);
        Vec3r p_end = Ts.back().block<3,1>(0,3);
        if (_draw_end_caps)
        {
            for (int i = 0; i < _cap_resolution; i++)
            {
                Real phi = (M_PI / 2.0) * i / _cap_resolution;

                Real ring_radius = robot->_r * std::cos(phi);
                Real z = robot->_r * std::sin(phi);

                for (int j = 0; j < _tubular_resolution; j++)
                {
                    Real theta = 2.0 * M_PI * j / _tubular_resolution;

                    Vec3r p(
                        ring_radius * std::cos(theta),
                        ring_radius * std::sin(theta),
                        z);

                    Vec3r transformed = R_end * p + p_end;
                    points->InsertNextPoint(transformed.data());
                }
            }
        }


        int total_num_rings = num_samples;
        if (_draw_end_caps)
            total_num_rings += 2*_cap_resolution;

        // create faces inbetween nodes
        vtkNew<vtkCellArray> faces;
        for (int si = 0; si < total_num_rings-1; si++)
        {
            for (int pi = 0; pi < _tubular_resolution; pi++)
            {
                // vertices for the quad face
                const int v1 = si*_tubular_resolution + pi;
                const int v2 = (pi != _tubular_resolution-1) ? v1 + 1 : si*_tubular_resolution;
                const int v3 = v2 + _tubular_resolution;
                const int v4 = v1 + _tubular_resolution;

                vtkNew<vtkTriangle> tri1;
                tri1->GetPointIds()->SetId(0, v1);
                tri1->GetPointIds()->SetId(1, v2);
                tri1->GetPointIds()->SetId(2, v3);

                vtkNew<vtkTriangle> tri2;
                tri2->GetPointIds()->SetId(0, v1);
                tri2->GetPointIds()->SetId(1, v3);
                tri2->GetPointIds()->SetId(2, v4);
                faces->InsertNextCell(tri1);
                faces->InsertNextCell(tri2);
            }
        }

        // end cap points
        if (_draw_end_caps)
        {
            Vec3r base_pt = R_beg * Vec3r(0, 0, -robot->_r) + p_beg;
            Vec3r end_pt = R_end * Vec3r(0, 0, robot->_r) + p_end;
            points->InsertNextPoint(base_pt.data());
            points->InsertNextPoint(end_pt.data());
        }
        else
        {
            points->InsertNextPoint( p_beg[0], p_beg[1], p_beg[2] );
            points->InsertNextPoint( p_end[0], p_end[1], p_end[2] );
        }

        // end cap faces
        // base faces
        for (int pi = 0; pi < _tubular_resolution; pi++)
        {
            const int v1 = pi;
            const int v2 = total_num_rings*_tubular_resolution;
            const int v3 = (pi != _tubular_resolution-1) ? v1 + 1 : 0;

            vtkNew<vtkTriangle> tri;
            tri->GetPointIds()->SetId(0, v1);
            tri->GetPointIds()->SetId(1, v2);
            tri->GetPointIds()->SetId(2, v3);
            faces->InsertNextCell(tri);
        }

        // tip faces
        for (int pi = 0; pi < _tubular_resolution; pi++)
        {
            const int v1 = (total_num_rings-1)*_tubular_resolution + pi;
            const int v2 = (pi != _tubular_resolution-1) ? v1 + 1 : (total_num_rings-1)*_tubular_resolution;
            const int v3 = total_num_rings*_tubular_resolution + 1;
            

            vtkNew<vtkTriangle> tri;
            tri->GetPointIds()->SetId(0, v1);
            tri->GetPointIds()->SetId(1, v2);
            tri->GetPointIds()->SetId(2, v3);
            faces->InsertNextCell(tri);
        }
        
        _vtk_poly_data->SetPoints(points);
        _vtk_poly_data->SetPolys(faces);
    }

    void _updatePolyData(const VecXr& q)
    {
        std::vector<Mat4r> Ts = robot->integratePoseAlongLength(_s_ips, q);

        vtkPoints* points = _vtk_poly_data->GetPoints();

        // base cap
        Mat3r R_beg = Ts.front().block<3,3>(0,0);
        Vec3r p_beg = Ts.front().block<3,1>(0,3);
        if (_draw_end_caps)
        {
            for (int i = _cap_resolution-1; i >= 0; i--)
            {
                Real phi = (M_PI / 2.0) * i / _cap_resolution;

                Real ring_radius = robot->_r * std::cos(phi);
                Real z = robot->_r * std::sin(phi);

                for (int j = 0; j < _tubular_resolution; j++)
                {
                    Real theta = 2.0 * M_PI * j / _tubular_resolution;

                    Vec3r p(
                        ring_radius * std::cos(theta),
                        ring_radius * std::sin(theta),
                        -z);

                    Vec3r transformed = R_beg * p + p_beg;
                    int offset = (_cap_resolution - 1 - i) * _tubular_resolution;
                    points->SetPoint(offset + j, transformed.data());
                }
            }
        }

        int cs_start_index = _draw_end_caps ? _cap_resolution*_tubular_resolution : 0;
        for (int si = 0; si < num_samples; si++)
        {
            Vec3r position = Ts[si].block<3,1>(0,3);
            Mat3r orientation = Ts[si].block<3,3>(0,0);
            // for (unsigned pi = 0; pi < cross_section_points.size(); pi++)
            for (int pi = 0; pi < _tubular_resolution; pi++)
            {
                // transform each point in local cross section to global frame
                const Vec3r transformed = orientation * _cross_section_points[pi] + position;
                points->SetPoint( cs_start_index + si*_tubular_resolution + pi, transformed.data() );
            }
        }

        // end cap
        Mat3r R_end = Ts.back().block<3,3>(0,0);
        Vec3r p_end = Ts.back().block<3,1>(0,3);
        if (_draw_end_caps)
        {
            for (int i = 0; i < _cap_resolution; i++)
            {
                Real phi = (M_PI / 2.0) * i / _cap_resolution;

                Real ring_radius = robot->_r * std::cos(phi);
                Real z = robot->_r * std::sin(phi);

                for (int j = 0; j < _tubular_resolution; j++)
                {
                    Real theta = 2.0 * M_PI * j / _tubular_resolution;

                    Vec3r p(
                        ring_radius * std::cos(theta),
                        ring_radius * std::sin(theta),
                        z);

                    Vec3r transformed = R_end * p + p_end;
                    points->SetPoint((_cap_resolution+num_samples)*_tubular_resolution + i*_tubular_resolution + j, transformed.data());
                }
            }
        }

        if (_draw_end_caps)
        {
            Vec3r base_pt = R_beg * Vec3r(0, 0, -robot->_r) + p_beg;
            Vec3r end_pt = R_end * Vec3r(0, 0, robot->_r) + p_end;
            points->SetPoint((2*_cap_resolution + num_samples)*_tubular_resolution, base_pt.data());
            points->SetPoint((2*_cap_resolution + num_samples)*_tubular_resolution+1, end_pt.data());
        }
        else
        {
            points->SetPoint(num_samples*_tubular_resolution, p_beg[0], p_beg[1], p_beg[2]);
            points->SetPoint(num_samples*_tubular_resolution+1, p_end[0], p_end[1], p_end[2]);
        }

        points->Modified();
    }

};


/////////////////////////////////////////////////////////////////////////////////////////////

/** Handles VTK window creation and rendering */
struct Visualization
{
    /** VTK renderer */
    vtkSmartPointer<vtkOpenGLRenderer> _renderer;
    /** VTK render window */
    vtkSmartPointer<vtkRenderWindow> _render_window;
    /** VTK render window interactor */
    vtkSmartPointer<vtkRenderWindowInteractor> _interactor;

    /** For saving frames as video */
    vtkSmartPointer<vtkWindowToImageFilter> _window_to_image;
    vtkSmartPointer<vtkPNGWriter> _png_writer;

    /** Graphics objects */
    std::vector<ContinuumRobotGraphicsObject> _graphics_objects;

    void displayWindow() { _render_window->Render(); }
    void interactorStart() { _interactor->Start(); }

    void addGraphicsObject(ContinuumRobot* robot, int num_samples)
    {
        _graphics_objects.emplace_back(robot, num_samples);
    }

    void setup()
    {
        // create renderer for actors in the scene
        _renderer = vtkSmartPointer<vtkOpenGLRenderer>::New();
        
        // add all the rods that may have been added before we set up
        for (auto& graphics_obj : _graphics_objects)
        {
            _renderer->AddActor(graphics_obj.actor());
        }
        

        _renderer->SetBackground(0.8, 0.8, 0.8);
        _renderer->SetAutomaticLightCreation(false);

        // set camera to look at (0,0,0) from its current position
        _renderer->ResetCamera();
        vtkCamera* camera = _renderer->GetActiveCamera();
        double focal_point[3];
        camera->GetFocalPoint(focal_point);

        // set focal point
        camera->SetFocalPoint(
            0,
            0,
            0
        );

        // position camera
        camera->SetPosition(
            1,
            1,
            1
        );

        camera->SetViewUp(0, 1, 0);
        camera->SetClippingRange(0.01, 1000.0);

        camera->SetParallelProjection(false);

        //////////////////////////////////////////////////////////
        // Create HDR lighting (if specified in the config)
        /////////////////////////////////////////////////////////

        // std::optional<std::string> hdr_filename = _render_config.hdrImageFilename();
        // if (hdr_filename.has_value())
        // {
        //     vtkNew<vtkTexture> hdr_texture;
        //     vtkNew<vtkHDRReader> reader;
        //     reader->SetFileName(hdr_filename.value().c_str());
        //     vtkNew<vtkImageShiftScale> scale;
        //     scale->SetInputConnection(reader->GetOutputPort());
        //     scale->SetScale(_render_config.hdrScaling()); // darker HDRI

        //     hdr_texture->SetInputConnection(scale->GetOutputPort());
        //     hdr_texture->SetColorModeToDirectScalars();
        //     hdr_texture->MipmapOn();
        //     hdr_texture->InterpolateOn();
        //     hdr_texture->SetInputConnection(scale->GetOutputPort());

        //     if (_render_config.createSkybox())
        //     {
        //         vtkNew<vtkSkybox> skybox;
        //         skybox->SetTexture(hdr_texture);
        //         skybox->SetFloorRight(0,0,1);
        //         skybox->SetProjection(vtkSkybox::Sphere);
        //         _renderer->AddActor(skybox);
        //     }

        //     _renderer->UseImageBasedLightingOn();
        //     _renderer->UseSphericalHarmonicsOff();
        //     _renderer->SetEnvironmentTexture(hdr_texture, false);
        // }
        

        ////////////////////////////////////////////////////////
        // Add lights
        ////////////////////////////////////////////////////////

        // vtkNew<vtkLight> light;
        // light->SetLightTypeToSceneLight();
        // light->SetPositional(true);
        // Vec3r light_pos = _render_config.sceneLightPosition();
        // Vec3r light_fp = _render_config.sceneLightFocalPoint();
        // light->SetPosition(light_pos[0], light_pos[1], light_pos[2]);
        // light->SetFocalPoint(light_fp[0], light_fp[1], light_fp[2]);
        // light->SetConeAngle(45);             
        // light->SetColor(1.0, 1.0, 1.0);
        // light->SetIntensity(1.0);
        // light->SetShadowAttenuation(0.5);
        // _renderer->AddLight(light);


        //////////////////////////////////////////////////////
        // Create the render window and interactor
        //////////////////////////////////////////////////////
        _render_window = vtkSmartPointer<vtkRenderWindow>::New();
        _render_window->AddRenderer(_renderer);
        _render_window->SetSize(800, 800);
        _render_window->SetWindowName("Simulation Playback");

        _interactor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
        vtkNew<vtkInteractorStyleTrackballCamera> style;
        // vtkNew<CustomVTKInteractorStyle> style;
        style->AutoAdjustCameraClippingRangeOff();
        _interactor->SetInteractorStyle(style);
        _interactor->SetRenderWindow(_render_window);
        
        
        
        /////////////////////////////////////////////////////////
        // Create the rendering passes and settings
        ////////////////////////////////////////////////////////
        _render_window->SetAlphaBitPlanes(1);
        // _render_window->SetMultiSamples(0);

        _renderer->SetUseDepthPeeling(true);
        _renderer->SetMaximumNumberOfPeels(50);
        _renderer->SetOcclusionRatio(0.1);

        // Core passes
        vtkNew<vtkRenderPassCollection> passes;

        vtkNew<vtkLightsPass> lightsPass;
        vtkNew<vtkOpaquePass> opaquePass;
        vtkNew<vtkTranslucentPass> translucentPass;

        // Shadow map
        vtkNew<vtkShadowMapPass> shadowPass;
        shadowPass->GetShadowMapBakerPass()->SetResolution(4096);

        // IMPORTANT: depth peeling must wrap translucent pass
        vtkNew<vtkDepthPeelingPass> peelingPass;
        peelingPass->SetTranslucentPass(translucentPass);
        peelingPass->SetMaximumNumberOfPeels(50);
        peelingPass->SetOcclusionRatio(0.1);

        // Order matters
        passes->AddItem(shadowPass->GetShadowMapBakerPass());
        passes->AddItem(lightsPass);
        passes->AddItem(opaquePass);
        passes->AddItem(shadowPass);
        passes->AddItem(peelingPass);

        vtkNew<vtkSequencePass> seq;
        seq->SetPasses(passes);

        vtkNew<vtkCameraPass> cameraP;
        cameraP->SetDelegatePass(seq);

        _renderer->SetPass(cameraP);

        // if rendering for video, set up the window-to-image and PNG writer
        // if (_render_config.renderForVideo())
        // {
        //     _window_to_image = vtkSmartPointer<vtkWindowToImageFilter>::New();
        //     _window_to_image->SetInput(_render_window);
        //     _window_to_image->SetInputBufferTypeToRGBA();
        //     _window_to_image->ReadFrontBufferOff();

        //     _png_writer = vtkSmartPointer<vtkPNGWriter>::New();
        //     _png_writer->SetInputConnection(_window_to_image->GetOutputPort());

        //     // enable offscreen rendering
        //     _render_window->SetOffScreenRendering(true);
        // }

        // vtkNew<vtkCallbackCommand> render_callback;
        // render_callback->SetCallback(GraphicsScene::renderCallback);
        // render_callback->SetClientData(this);
        _interactor->Initialize();
        // _interactor->AddObserver(vtkCommand::TimerEvent, render_callback);
        // _interactor->CreateRepeatingTimer(5);
    }
};