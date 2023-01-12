#include "ViewerApplication.hpp"

#include <iostream>
#include <numeric>

#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/quaternion.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtx/io.hpp>

#include "utils/cameras.hpp"

#include <stb_image_write.h>


void keyCallback(
    GLFWwindow *window, int key, int scancode, int action, int mods)
{
  if (key == GLFW_KEY_ESCAPE && action == GLFW_RELEASE) {
    glfwSetWindowShouldClose(window, 1);
  }
}

//////////////////////////////// Loading the glTF file ////////////////////////////////////
bool ViewerApplication::loadGltfFile(tinygltf::Model & model) 
{
  tinygltf::TinyGLTF loader;
  std::string err;
  std::string warn;

  bool ret = loader.LoadASCIIFromFile(&model, &err, &warn, m_gltfFilePath.string());

  if (!warn.empty()) {
    printf("Warn: %s\n", warn.c_str());
  }

  if (!err.empty()) {
    printf("Err: %s\n", err.c_str());
  }

  if (!ret) {
    printf("Failed to parse glTF\n");
    return false;
  }

  return ret;
}

//////////////////////////////// Creation of Buffer Objects ////////////////////////////////////
std::vector<GLuint> ViewerApplication::createBufferObjects( const tinygltf::Model &model)
{
  std::vector<GLuint> bufferObjects(model.buffers.size(), 0); 
  
  glGenBuffers(model.buffers.size(), bufferObjects.data());
  for (size_t i = 0; i < model.buffers.size(); ++i) {
    glBindBuffer(GL_ARRAY_BUFFER, bufferObjects[i]);
    glBufferStorage(GL_ARRAY_BUFFER, model.buffers[i].data.size(), model.buffers[i].data.data(), 0);
  }
  glBindBuffer(GL_ARRAY_BUFFER, 0);
  return bufferObjects;
}


//////////////////////////////// Creation of Vertex Array Objects ////////////////////////////////////
std::vector<GLuint> ViewerApplication::createVertexArrayObjects( const tinygltf::Model &model, const std::vector<GLuint> &bufferObjects, std::vector<VaoRange> & meshIndexToVaoRange) 
{
  std::vector<GLuint> vertexArrayObjects;

  const GLuint VERTEX_ATTRIB_POSITION_IDX = 0;
  const GLuint VERTEX_ATTRIB_NORMAL_IDX = 1;
  const GLuint VERTEX_ATTRIB_TEXCOORD0_IDX = 2;

  meshIndexToVaoRange.resize(model.meshes.size());

  for(int meshIdx = 0; meshIdx < model.meshes.size(); meshIdx++) 
  { 
    const auto &mesh = model.meshes[meshIdx];

    auto &vaoRange = meshIndexToVaoRange[meshIdx];
    vaoRange.begin = GLsizei(vertexArrayObjects.size()); 
    vaoRange.count = GLsizei(mesh.primitives.size()); 

    vertexArrayObjects.resize(vertexArrayObjects.size() + mesh.primitives.size());

    glGenVertexArrays(vaoRange.count, &vertexArrayObjects[vaoRange.begin]);  
    for (auto primitiveIdx = 0; primitiveIdx < mesh.primitives.size(); primitiveIdx++)
    {
      const auto vao = vertexArrayObjects[vaoRange.begin + primitiveIdx];
      const auto &primitive = mesh.primitives[primitiveIdx];

      glBindVertexArray(vao);

      { // I'm opening a scope because I want to reuse the variable iterator in the code for NORMAL and TEXCOORD_0
        const auto iterator = primitive.attributes.find("POSITION");
        if (iterator != end(primitive.attributes)) { // If "POSITION" has been found in the map
          // (*iterator).first is the key "POSITION", (*iterator).second is the value, ie. the index of the accessor for this attribute
          const auto accessorIdx = (*iterator).second;
          const auto &accessor = model.accessors[accessorIdx];              // TODO get the correct tinygltf::Accessor from model.accessors
          const auto &bufferView = model.bufferViews[accessor.bufferView]; // TODO get the correct tinygltf::BufferView from model.bufferViews. You need to use the accessor
          const auto bufferIdx = bufferView.buffer;                       // TODO get the index of the buffer used by the bufferView (you need to use it)
          const auto bufferObject = bufferObjects[bufferIdx];            // TODO get the correct buffer object from the buffer index

          // TODO Enable the vertex attrib array corresponding to POSITION with glEnableVertexAttribArray (you need to use VERTEX_ATTRIB_POSITION_IDX which has to be defined at the top of the cpp file)
          glEnableVertexAttribArray(VERTEX_ATTRIB_POSITION_IDX);
          // TODO Bind the buffer object to GL_ARRAY_BUFFER
          glBindBuffer(GL_ARRAY_BUFFER, bufferObject);

          const auto byteOffset = accessor.byteOffset + bufferView.byteOffset; // TODO Compute the total byte offset using the accessor and the buffer view
          // TODO Call glVertexAttribPointer with the correct arguments.
          // Remember size is obtained with accessor.type, type is obtained with accessor.componentType.
          // The stride is obtained in the bufferView, normalized is always GL_FALSE, and pointer is the byteOffset (don't forget the cast).
          glVertexAttribPointer(VERTEX_ATTRIB_POSITION_IDX, accessor.type, accessor.componentType, 
                                GL_FALSE, GLsizei(bufferView.byteStride), (const GLvoid *)byteOffset);
        };
      }

      { 
        const auto iterator = primitive.attributes.find("NORMAL");
        if (iterator != end(primitive.attributes)) { 
          const auto accessorIdx = (*iterator).second;
          const auto &accessor = model.accessors[accessorIdx];   
          const auto &bufferView = model.bufferViews[accessor.bufferView]; 
          const auto bufferIdx = bufferView.buffer;                       
          const auto bufferObject = bufferObjects[bufferIdx];            

          glEnableVertexAttribArray(VERTEX_ATTRIB_NORMAL_IDX);
          glBindBuffer(GL_ARRAY_BUFFER, bufferObject);

          const auto byteOffset = accessor.byteOffset + bufferView.byteOffset;
          glVertexAttribPointer(VERTEX_ATTRIB_NORMAL_IDX, accessor.type, accessor.componentType, 
                                GL_FALSE, GLsizei(bufferView.byteStride), (const GLvoid *)byteOffset);
        };
      }

      { 
        const auto iterator = primitive.attributes.find("TEXCOORD_0");
        if (iterator != end(primitive.attributes)) { 
          const auto accessorIdx = (*iterator).second;
          const auto &accessor = model.accessors[accessorIdx];   
          const auto &bufferView = model.bufferViews[accessor.bufferView]; 
          const auto bufferIdx = bufferView.buffer;                       
          const auto bufferObject = bufferObjects[bufferIdx];            

          glEnableVertexAttribArray(VERTEX_ATTRIB_TEXCOORD0_IDX);
          glBindBuffer(GL_ARRAY_BUFFER, bufferObject);

          const auto byteOffset = accessor.byteOffset + bufferView.byteOffset;
          glVertexAttribPointer(VERTEX_ATTRIB_TEXCOORD0_IDX, accessor.type, accessor.componentType, 
                                GL_FALSE, GLsizei(bufferView.byteStride), (const GLvoid *)byteOffset);
        };
      }

      if (primitive.indices >= 0) {
        const auto accessorIdx = primitive.indices;
        const auto &accessor = model.accessors[accessorIdx];
        const auto &bufferView = model.bufferViews[accessor.bufferView];
        const auto bufferIdx = bufferView.buffer;
        const auto bufferObject = bufferObjects[bufferIdx];

        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, bufferObject);
      }
    }
  }

  return vertexArrayObjects;
}

int ViewerApplication::run()
{
  // Loader shaders
  const auto glslProgram =
      compileProgram({m_ShadersRootPath / m_vertexShader,
          m_ShadersRootPath / m_fragmentShader});

  const auto modelViewProjMatrixLocation = glGetUniformLocation(glslProgram.glId(), "uModelViewProjMatrix");
  const auto modelViewMatrixLocation     = glGetUniformLocation(glslProgram.glId(), "uModelViewMatrix");
  const auto normalMatrixLocation        = glGetUniformLocation(glslProgram.glId(), "uNormalMatrix");
  
  tinygltf::Model model;
  if (!loadGltfFile(model)) {
    return -1;
  }
  glm::vec3 bboxMin, bboxMax;
  computeSceneBounds(model, bboxMin, bboxMax);

  const auto diag = bboxMax - bboxMin;
  auto maxDistance = glm::length(diag);
  const auto projMatrix =
      glm::perspective(70.f, float(m_nWindowWidth) / m_nWindowHeight,
          0.001f * maxDistance, 1.5f * maxDistance);

  // TODO Implement a new CameraController model and use it instead. Propose the choice from the GUI
  std::unique_ptr<CameraController> cameraController =  std::make_unique<TrackballCameraController>(m_GLFWHandle.window(), 0.5f * maxDistance);
  if (m_hasUserCamera) {
    cameraController->setCamera(m_userCamera);
  } else {
    const auto center = 0.5f * (bboxMax + bboxMin);
    const auto up = glm::vec3(0, 1, 0);
    const auto eye = diag.z > 0 ? center + diag : center + 2.f * glm::cross(diag, up);
    cameraController->setCamera(Camera{eye, center, up});
  }

  // Loading the glTF file
  loadGltfFile(model);

  // Creation of Buffer Objects
  const auto bufferObjects = createBufferObjects(model);

  // Creation of Vertex Array Objects
  std::vector<VaoRange> meshToVertexArrays;
  const auto vertexArrayObjects = createVertexArrayObjects(model, bufferObjects, meshToVertexArrays);

  // Setup OpenGL state for rendering
  glEnable(GL_DEPTH_TEST);
  glslProgram.use();

  // Lambda function to draw the scene
  const auto drawScene = [&](const Camera &camera) {
    glViewport(0, 0, m_nWindowWidth, m_nWindowHeight);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    const auto viewMatrix = camera.getViewMatrix();

    // The recursive function that should draw a node
    // We use a std::function because a simple lambda cannot be recursive
    const std::function<void(int, const glm::mat4 &)> drawNode =
        [&](int nodeIdx, const glm::mat4 &parentMatrix) {
          // The drawNode function
          const auto &node = model.nodes[nodeIdx];
          glm::mat4 modelMatrix = getLocalToWorldMatrix(node, parentMatrix);

          if(node.mesh >= 0)
          {
            const auto mvMatrix  = viewMatrix * modelMatrix;
            const auto mvpMatrix = projMatrix * mvMatrix;

            glUniformMatrix4fv(modelViewProjMatrixLocation, 1, GL_FALSE, glm::value_ptr(mvpMatrix));
            glUniformMatrix4fv(modelViewMatrixLocation,     1, GL_FALSE, glm::value_ptr(mvMatrix));         
            glUniformMatrix4fv(normalMatrixLocation,        1, GL_FALSE, glm::value_ptr(glm::transpose(glm::inverse(mvMatrix))));
            
            const auto &mesh = model.meshes[node.mesh];
            const auto &vaoRange = meshToVertexArrays[node.mesh];



            for (auto primIdx = 0; primIdx < mesh.primitives.size(); primIdx++) {
              const auto vao = vertexArrayObjects[vaoRange.begin + primIdx];
              const auto &primitive = mesh.primitives[primIdx];

              glBindVertexArray(vao);

              if (primitive.indices >= 0) 
              {
                const auto &accessor   = model.accessors[primitive.indices];
                const auto &bufferView = model.bufferViews[accessor.bufferView];
                const auto byteOffset  = accessor.byteOffset + bufferView.byteOffset;
                
                glDrawElements(primitive.mode, GLsizei(accessor.count), accessor.componentType, (const GLvoid *) byteOffset);
              } 
              
              else 
              {
                const auto accessorIdx = (*begin(primitive.attributes)).second;
                const auto &accessor = model.accessors[accessorIdx];
                glDrawArrays(primitive.mode, 0, GLsizei(accessor.count));
              }
            }
          }

          // Draw children
          for (const auto childNodeIdx : node.children) {
            drawNode(childNodeIdx, modelMatrix);
          }
        };

    // Draw the scene referenced by gltf file
    if (model.defaultScene >= 0) {
      // Draw all nodes
      for(auto node : model.scenes[model.defaultScene].nodes)
      {
        drawNode(node, glm::mat4(1));
      }
    }
  };

  if(!m_OutputPath.empty())
  {
    const auto numComponents = 3;

    std::vector<unsigned char> pixels(m_nWindowWidth * m_nWindowHeight * numComponents);

    renderToImage(m_nWindowWidth, m_nWindowHeight, numComponents, pixels.data(), [&]() {
      drawScene(cameraController->getCamera());
    });

    flipImageYAxis(
        m_nWindowWidth, m_nWindowHeight, numComponents, pixels.data());

    // Write png on disk
    const auto strPath = m_OutputPath.string();
    stbi_write_png(
        strPath.c_str(), m_nWindowWidth, m_nWindowHeight, numComponents, pixels.data(), 0);

    return 0; 
  }

  // Loop until the user closes the window
  for (auto iterationCount = 0u; !m_GLFWHandle.shouldClose();
       ++iterationCount) {
    const auto seconds = glfwGetTime();

    const auto camera = cameraController->getCamera();
    drawScene(camera);

    // GUI code:
    imguiNewFrame();

    {
      ImGui::Begin("GUI");
      ImGui::Text("Application average %.3f ms/frame (%.1f FPS)",
          1000.0f / ImGui::GetIO().Framerate, ImGui::GetIO().Framerate);
      if (ImGui::CollapsingHeader("Camera", ImGuiTreeNodeFlags_DefaultOpen)) {
        ImGui::Text("eye: %.3f %.3f %.3f", camera.eye().x, camera.eye().y,
            camera.eye().z);
        ImGui::Text("center: %.3f %.3f %.3f", camera.center().x,
            camera.center().y, camera.center().z);
        ImGui::Text(
            "up: %.3f %.3f %.3f", camera.up().x, camera.up().y, camera.up().z);

        ImGui::Text("front: %.3f %.3f %.3f", camera.front().x, camera.front().y,
            camera.front().z);
        ImGui::Text("left: %.3f %.3f %.3f", camera.left().x, camera.left().y,
            camera.left().z);

        if (ImGui::Button("CLI camera args to clipboard")) {
          std::stringstream ss;
          ss << "--lookat " << camera.eye().x << "," << camera.eye().y << ","
             << camera.eye().z << "," << camera.center().x << ","
             << camera.center().y << "," << camera.center().z << ","
             << camera.up().x << "," << camera.up().y << "," << camera.up().z;
          const auto str = ss.str();
          glfwSetClipboardString(m_GLFWHandle.window(), str.c_str());
        }
        
        static int cameraControllerType = 0;
        const auto cameraControllerTypeChanged =
            ImGui::RadioButton("Trackball", &cameraControllerType, 0) ||
            ImGui::RadioButton("First Person", &cameraControllerType, 1);
        if (cameraControllerTypeChanged) {
          const auto currentCamera = cameraController->getCamera();
          if (cameraControllerType == 0) {
            cameraController = std::make_unique<TrackballCameraController>(
                m_GLFWHandle.window(), 0.5f * maxDistance);
          } else {
            cameraController = std::make_unique<FirstPersonCameraController>(
                m_GLFWHandle.window(), 0.5f * maxDistance);
          }
          cameraController->setCamera(currentCamera);
        }

      }
      ImGui::End();
    }

    imguiRenderFrame();

    glfwPollEvents(); // Poll for and process events

    auto ellapsedTime = glfwGetTime() - seconds;
    auto guiHasFocus =
        ImGui::GetIO().WantCaptureMouse || ImGui::GetIO().WantCaptureKeyboard;
    if (!guiHasFocus) {
      cameraController->update(float(ellapsedTime));
    }

    m_GLFWHandle.swapBuffers(); // Swap front and back buffers
  }

  // TODO clean up allocated GL data

  return 0;
}

ViewerApplication::ViewerApplication(const fs::path &appPath, uint32_t width,
    uint32_t height, const fs::path &gltfFile,
    const std::vector<float> &lookatArgs, const std::string &vertexShader,
    const std::string &fragmentShader, const fs::path &output) :
    m_nWindowWidth(width),
    m_nWindowHeight(height),
    m_AppPath{appPath},
    m_AppName{m_AppPath.stem().string()},
    m_ImGuiIniFilename{m_AppName + ".imgui.ini"},
    m_ShadersRootPath{m_AppPath.parent_path() / "shaders"},
    m_gltfFilePath{gltfFile},
    m_OutputPath{output}
{
  if (!lookatArgs.empty()) {
    m_hasUserCamera = true;
    m_userCamera =
        Camera{glm::vec3(lookatArgs[0], lookatArgs[1], lookatArgs[2]),
            glm::vec3(lookatArgs[3], lookatArgs[4], lookatArgs[5]),
            glm::vec3(lookatArgs[6], lookatArgs[7], lookatArgs[8])};
  }

  if (!vertexShader.empty()) {
    m_vertexShader = vertexShader;
  }

  if (!fragmentShader.empty()) {
    m_fragmentShader = fragmentShader;
  }

  ImGui::GetIO().IniFilename =
      m_ImGuiIniFilename.c_str(); // At exit, ImGUI will store its windows
                                  // positions in this file

  glfwSetKeyCallback(m_GLFWHandle.window(), keyCallback);

  printGLVersion();
}
