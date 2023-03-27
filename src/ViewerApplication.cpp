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
std::vector<GLuint> ViewerApplication::createVertexArrayObjects( const tinygltf::Model &model, const std::vector<GLuint> &bufferObjects, const GLuint bufferTangentObject, std::vector<VaoRange> & meshIndexToVaoRange) 
{
  std::vector<GLuint> vertexArrayObjects;

  const GLuint VERTEX_ATTRIB_POSITION_IDX  = 0;
  const GLuint VERTEX_ATTRIB_NORMAL_IDX    = 1;
  const GLuint VERTEX_ATTRIB_TEXCOORD0_IDX = 2;
  const GLuint VERTEX_ATTRIB_TANGENT_IDX   = 3; 
  const GLuint VERTEX_ATTRIB_BITANGENT_IDX = 4; 
  
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
          const auto &accessor = model.accessors[accessorIdx];              // get the correct tinygltf::Accessor from model.accessors
          const auto &bufferView = model.bufferViews[accessor.bufferView]; // get the correct tinygltf::BufferView from model.bufferViews. You need to use the accessor
          const auto bufferIdx = bufferView.buffer;                       // get the index of the buffer used by the bufferView (you need to use it)
          const auto bufferObject = bufferObjects[bufferIdx];            // get the correct buffer object from the buffer index

          // Enable the vertex attrib array corresponding to POSITION with glEnableVertexAttribArray (you need to use VERTEX_ATTRIB_POSITION_IDX which has to be defined at the top of the cpp file)
          glEnableVertexAttribArray(VERTEX_ATTRIB_POSITION_IDX);
          // Bind the buffer object to GL_ARRAY_BUFFER
          glBindBuffer(GL_ARRAY_BUFFER, bufferObject);

          const auto byteOffset = accessor.byteOffset + bufferView.byteOffset; // Compute the total byte offset using t0n't forget the cast).
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
      
      // TANGENTS 
      glEnableVertexAttribArray(VERTEX_ATTRIB_TANGENT_IDX);
      glBindBuffer(GL_ARRAY_BUFFER, bufferTangentObject);
      glVertexAttribPointer(VERTEX_ATTRIB_TANGENT_IDX, 3, GL_FLOAT, GL_FALSE, 0, (const GLvoid *)0);
      
      // BITANGENTS
      /*
      glEnableVertexAttribArray(VERTEX_ATTRIB_BITANGENT_IDX);
      glBindBuffer(GL_ARRAY_BUFFER, bufferBitangentObject);
      glVertexAttribPointer(VERTEX_ATTRIB_BITANGENT_IDX, 3, GL_FLOAT, GL_FALSE, 0, (const GLvoid *)0);
      */

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

glm::vec3 ViewerApplication::createTangentBitangent(const glm::vec3 * position, const glm::vec2 * texture)
{ 
  glm::vec3 tangent;

  glm::vec3 edge1 = position[1] - position[0];
  glm::vec3 edge2 = position[2] - position[0];
  glm::vec2 deltaUV1 = texture[1] - texture[0];
  glm::vec2 deltaUV2 = texture[2] - texture[0];  

  float f = 1.0f / (deltaUV1.x * deltaUV2.y - deltaUV2.x * deltaUV1.y);

  tangent.x = f * (deltaUV2.y * edge1.x - deltaUV1.y * edge2.x);
  tangent.y = f * (deltaUV2.y * edge1.y - deltaUV1.y * edge2.y);
  tangent.z = f * (deltaUV2.y * edge1.z - deltaUV1.y * edge2.z);

  return tangent;
/*
  bitangent.x = f * (-deltaUV2.x * edge1.x + deltaUV1.x * edge2.x);
  bitangent.y = f * (-deltaUV2.x * edge1.y + deltaUV1.x * edge2.y);
  bitangent.z = f * (-deltaUV2.x * edge1.z + deltaUV1.x * edge2.z);
  */
}

//////////////////////////////// Creation of Tangent Buffer Objects ////////////////////////////////////
GLuint createBufferFromArrayObjects(glm::vec3 * array, int arraySize)
{
  GLuint bufferObject; 
  glGenBuffers((GLsizei) 1, &bufferObject);
  glBindBuffer(GL_ARRAY_BUFFER, bufferObject);
  glBufferStorage(GL_ARRAY_BUFFER, arraySize * sizeof(glm::vec3), array, GL_STATIC_DRAW);
  glBindBuffer(GL_ARRAY_BUFFER, 0);
  return bufferObject;
}

GLuint ViewerApplication::createTangentBitangentArrays(const tinygltf::Model &model)
{
  glm::vec3 * tangents = nullptr;
  int tangentNumber;

  if (model.defaultScene >= 0) {
    const std::function<void(int, const glm::mat4 &)> updateBounds =
        [&](int nodeIdx, const glm::mat4 &parentMatrix) {
          const auto &node = model.nodes[nodeIdx];
          const glm::mat4 modelMatrix = getLocalToWorldMatrix(node, parentMatrix);

          glm::vec3 localPositions[3];
          glm::vec2 localTextures[3];
          glm::vec3 tangent;
          //glm::vec3 bitangent;

          if (node.mesh >= 0) {
            const auto &mesh = model.meshes[node.mesh];
            
            // TODO, Passer les tableau en tableau de tableau au cas ou mesh.primites est sup a 1
            //tangents = (glm::vec3*) malloc(model.buffers[0].data.size() * sizeof(glm::vec3));
            //bitangents = (glm::vec3*) malloc(model.buffers[0].data.size() * sizeof(glm::vec3));
            //tangents.resize(model.buffers[0].data.size());

            for (size_t pIdx = 0; pIdx < mesh.primitives.size(); ++pIdx) {
              const auto &primitive        = mesh.primitives[pIdx];             
              const auto positionAttrIdxIt = primitive.attributes.find("POSITION");
              const auto textureAttrIdxIt  = primitive.attributes.find("TEXCOORD_0");

              if (positionAttrIdxIt == end(primitive.attributes) || textureAttrIdxIt == end(primitive.attributes)) {
                continue;
              }

              // Accessors
              const auto &positionAccessor = model.accessors[(*positionAttrIdxIt).second];
              const auto &textureAccessor  = model.accessors[(*textureAttrIdxIt).second];

              if (positionAccessor.type != 3) {
                std::cerr << "Position accessor with type != VEC3, skipping" << std::endl;
                continue;
              } if (textureAccessor.type != 2) {
                std::cerr << "Texture accessor with type != VEC2, skipping" << std::endl;
                continue;
              }

              // Position informations
              const auto &positionBufferView = model.bufferViews[positionAccessor.bufferView];
              const auto positionByteOffset  = positionAccessor.byteOffset + positionBufferView.byteOffset;
              const auto &positionBuffer     = model.buffers[positionBufferView.buffer];
              const auto positionByteStride  = positionBufferView.byteStride ? positionBufferView.byteStride : 3 * sizeof(float);
              // Texture informations
              const auto &textureBufferView = model.bufferViews[textureAccessor.bufferView];
              const auto textureByteOffset  = textureAccessor.byteOffset + textureBufferView.byteOffset;
              const auto &textureBuffer     = model.buffers[textureBufferView.buffer];
              const auto textureByteStride  = textureBufferView.byteStride ? textureBufferView.byteStride : 2 * sizeof(float);

              if(tangents == nullptr) {
                tangents = new glm::vec3[positionBuffer.data.size()];
              }

              if (primitive.indices >= 0) {
                const auto &indexAccessor   = model.accessors[primitive.indices];
                const auto &indexBufferView = model.bufferViews[indexAccessor.bufferView];
                const auto indexByteOffset  = indexAccessor.byteOffset + indexBufferView.byteOffset;
                const auto &indexBuffer     = model.buffers[indexBufferView.buffer];
                auto indexByteStride        = indexBufferView.byteStride;

                switch (indexAccessor.componentType) {
                  default:
                    std::cerr << "Primitive index accessor with bad componentType " << indexAccessor.componentType << ", skipping it." << std::endl;
                    continue;
                  case TINYGLTF_COMPONENT_TYPE_UNSIGNED_BYTE:
                    indexByteStride = indexByteStride ? indexByteStride : sizeof(uint8_t);
                    break;
                  case TINYGLTF_COMPONENT_TYPE_UNSIGNED_SHORT:
                    indexByteStride = indexByteStride ? indexByteStride : sizeof(uint16_t);
                    break;
                  case TINYGLTF_COMPONENT_TYPE_UNSIGNED_INT:
                    indexByteStride = indexByteStride ? indexByteStride : sizeof(uint32_t);
                    break;
                }

                uint32_t index[3];
                size_t i, j;

                for (i = 0; i < indexAccessor.count; i += 3) {        
                  for (j = 0; j < 3; j++) {
                    switch (indexAccessor.componentType) {
                      case TINYGLTF_COMPONENT_TYPE_UNSIGNED_BYTE:
                        index[j] = *((const uint8_t  *) &indexBuffer.data[indexByteOffset + indexByteStride * (i + j)]);
                        break;
                      case TINYGLTF_COMPONENT_TYPE_UNSIGNED_SHORT:
                        index[j] = *((const uint16_t *) &indexBuffer.data[indexByteOffset + indexByteStride * (i + j)]);
                        break;
                      case TINYGLTF_COMPONENT_TYPE_UNSIGNED_INT:
                        index[j] = *((const uint32_t *) &indexBuffer.data[indexByteOffset + indexByteStride * (i + j)]);
                        break;
                    }
                    localPositions[j] = *((const glm::vec3 *) &positionBuffer.data[positionByteOffset + positionByteStride * index[j]]);
                    localTextures[j] = *((const glm::vec2 *) &textureBuffer.data[textureByteOffset + textureByteStride * index[j]]);
                  }

                  tangent = createTangentBitangent(localPositions, localTextures);

                  ///////////////
                  for (j = 0; j < 3; j++) {
                    tangentNumber++;
                    tangents[positionByteOffset + positionByteStride * index[j]] = tangent;
                  }
                }
              } else {
                size_t i, j;
                for (size_t i = 0; i < positionAccessor.count; i += 3) {
                  for (size_t j = 0; j < 3; j++) {
                    localPositions[j] = *((const glm::vec3 *) &positionBuffer.data[positionByteOffset + positionByteStride * (i + j)]);
                    localTextures[j] = *((const glm::vec2 *) &textureBuffer.data[textureByteOffset + textureByteStride * (i + j)]);
                  }

                  tangent = createTangentBitangent(localPositions, localTextures);

                  ///////////////
                  for (j = 0; j < 3; j++) {
                    tangents[positionByteOffset + positionByteStride * (i + j)] = tangent;
                    tangentNumber++;
                  }
                }
              }
            }
          }
          
          for (const auto childNodeIdx : node.children) {
            updateBounds(childNodeIdx, modelMatrix);
          }
          
        };

    
    for (const auto nodeIdx : model.scenes[model.defaultScene].nodes) {
      updateBounds(nodeIdx, glm::mat4(1));
    }
  }

  return createBufferFromArrayObjects(tangents, tangentNumber);
}


std::vector<GLuint> ViewerApplication::createTextureObjects(const tinygltf::Model &model) const 
{
  std::vector<GLuint> textureObjects(model.textures.size(), 0);

  tinygltf::Sampler defaultSampler;
  defaultSampler.minFilter = GL_LINEAR;
  defaultSampler.magFilter = GL_LINEAR;
  defaultSampler.wrapS = GL_REPEAT;
  defaultSampler.wrapT = GL_REPEAT;
  defaultSampler.wrapR = GL_REPEAT;

  glActiveTexture(GL_TEXTURE0);

  glGenTextures(GLsizei(model.textures.size()), textureObjects.data());
  for (size_t i = 0; i < model.textures.size(); ++i) {
    const auto &texture = model.textures[i];
    assert(texture.source >= 0);
    const auto &image = model.images[texture.source];

    const auto &sampler = texture.sampler >= 0 ? model.samplers[texture.sampler] : defaultSampler;
    glBindTexture(GL_TEXTURE_2D, textureObjects[i]);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, image.width, image.height, 0, GL_RGBA, image.pixel_type, image.image.data());
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, sampler.minFilter != -1 ? sampler.minFilter : GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, sampler.magFilter != -1 ? sampler.magFilter : GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, sampler.wrapS);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, sampler.wrapT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_R, sampler.wrapR);

    if (sampler.minFilter == GL_NEAREST_MIPMAP_NEAREST || sampler.minFilter == GL_NEAREST_MIPMAP_LINEAR ||
        sampler.minFilter == GL_LINEAR_MIPMAP_NEAREST  || sampler.minFilter == GL_LINEAR_MIPMAP_LINEAR) {
      glGenerateMipmap(GL_TEXTURE_2D);
    }
  }
  glBindTexture(GL_TEXTURE_2D, 0);

  return textureObjects;
}

int ViewerApplication::run()
{
  // Loader shaders
  const auto glslProgram = compileProgram({m_ShadersRootPath / m_vertexShader, m_ShadersRootPath / m_fragmentShader});

  const auto modelViewProjMatrixLocation       = glGetUniformLocation(glslProgram.glId(), "uModelViewProjMatrix");
  const auto modelViewMatrixLocation           = glGetUniformLocation(glslProgram.glId(), "uModelViewMatrix");
  const auto normalMatrixLocation              = glGetUniformLocation(glslProgram.glId(), "uNormalMatrix");
  const auto modelMatrixLocation               = glGetUniformLocation(glslProgram.glId(), "uModelMatrix");

  const auto lightDirectionLocation            = glGetUniformLocation(glslProgram.glId(), "uLightDirection");
  const auto lightIntensityLocation            = glGetUniformLocation(glslProgram.glId(), "uLightIntensity");

  const auto baseColorTextureLocation          = glGetUniformLocation(glslProgram.glId(), "uBaseColorTexture");
  const auto baseColorFactorLocation           = glGetUniformLocation(glslProgram.glId(), "uBaseColorFactor");

  const auto metallicRoughnessTextureLocation  = glGetUniformLocation(glslProgram.glId(), "uMetallicRoughnessTexture");
  const auto metallicFactorLocation            = glGetUniformLocation(glslProgram.glId(), "uMetallicFactor");
  const auto roughnessFactorLocation           = glGetUniformLocation(glslProgram.glId(), "uRoughnessFactor");

  const auto emissiveTextureLocation           = glGetUniformLocation(glslProgram.glId(), "uEmissiveTexture");
  const auto emissiveFactorLocation            = glGetUniformLocation(glslProgram.glId(), "uEmissiveFactor");

  const auto occlusionTextureLocation          = glGetUniformLocation(glslProgram.glId(), "uOcclusionTexture");
  const auto occlusionStrengthLocation         = glGetUniformLocation(glslProgram.glId(), "uOcclusionStrength");
  const auto applyOcclusionLocation            = glGetUniformLocation(glslProgram.glId(), "uApplyOcclusion");
  const auto normalMapLocation                 = glGetUniformLocation(glslProgram.glId(), "uNormalMap");
  
  const auto viewNormalMapLocation             = glGetUniformLocation(glslProgram.glId(), "uViewNormalMap");
  
  tinygltf::Model model;
  if (!loadGltfFile(model)) {
    return -1;
  }

  // Set light parameters
  glm::vec3 lightDirection(1, 1, 1);
  glm::vec3 lightIntensity(1, 1, 1);
  bool lightFromCamera = false;
  bool applyOcclusion = true;
  bool viewNormalMap = false;

  // Creation of Texture Objects
  const auto textureObjects = createTextureObjects(model);
  GLuint whiteTexture = 0;

  // Create white texture for object with no base color texture
  glGenTextures(1, &whiteTexture);
  glBindTexture(GL_TEXTURE_2D, whiteTexture);
  float white[] = {1, 1, 1, 1};
  glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, 1, 1, 0, GL_RGBA, GL_FLOAT, white);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_R, GL_REPEAT);
  glBindTexture(GL_TEXTURE_2D, 0);

  ///////
  glm::vec3 bboxMin, bboxMax;
  computeSceneBounds(model, bboxMin, bboxMax);
  
  const auto diag = bboxMax - bboxMin;
  auto maxDistance = glm::length(diag);
  const auto projMatrix = glm::perspective(70.f, float(m_nWindowWidth) / m_nWindowHeight, 0.001f * maxDistance, 1.5f * maxDistance);

  // Implement a new CameraController model and use it instead. Propose the choice from the GUI
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

  // Creation of two Vertex Array of tangent and bitangent 
  const auto bufferTangentsObject = createTangentBitangentArrays(model); 

  // Creation of Vertex Array Objects
  std::vector<VaoRange> meshToVertexArrays;
  const auto vertexArrayObjects = createVertexArrayObjects(model, bufferObjects, bufferTangentsObject, meshToVertexArrays);

  // Setup OpenGL state for rendering
  glEnable(GL_DEPTH_TEST);
  glslProgram.use();

  const auto bindMaterial = [&](const auto materialIndex) {
    if (materialIndex >= 0) {
      const auto &material = model.materials[materialIndex];
      const auto &pbrMetallicRoughness = material.pbrMetallicRoughness;

      //// Base Color Factor ////
      if (baseColorFactorLocation >= 0) {
        glUniform4f(baseColorFactorLocation,
            (float)pbrMetallicRoughness.baseColorFactor[0],
            (float)pbrMetallicRoughness.baseColorFactor[1],
            (float)pbrMetallicRoughness.baseColorFactor[2],
            (float)pbrMetallicRoughness.baseColorFactor[3]);
      }

      //// Base Color Texture ////
      if (baseColorTextureLocation >= 0) {
        auto textureObject = whiteTexture;
        if (pbrMetallicRoughness.baseColorTexture.index >= 0) {
          const auto &texture =
              model.textures[pbrMetallicRoughness.baseColorTexture.index];
          if (texture.source >= 0) {
            textureObject = textureObjects[texture.source];
          }
        }
        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_2D, textureObject);
        glUniform1i(baseColorTextureLocation, 0);
      }

      //// Metallic Factor ////
      if (metallicFactorLocation >= 0) {
        glUniform1f(metallicFactorLocation, (float)pbrMetallicRoughness.metallicFactor);
      }

      //// Roughness Factor ////
      if (roughnessFactorLocation >= 0) {
        glUniform1f(roughnessFactorLocation, (float)pbrMetallicRoughness.roughnessFactor);
      }

      //// Metallic Roughness Texture ////
      if (metallicRoughnessTextureLocation >= 0) {
        auto textureObject = 0u;
        if (pbrMetallicRoughness.metallicRoughnessTexture.index >= 0) {
          const auto &texture = model.textures[pbrMetallicRoughness.metallicRoughnessTexture.index];
          if (texture.source >= 0) {
            textureObject = textureObjects[texture.source];
          }
        }
        glActiveTexture(GL_TEXTURE1);
        glBindTexture(GL_TEXTURE_2D, textureObject);
        glUniform1i(metallicRoughnessTextureLocation, 1);
      }

      //// Emissive Factor ////
      if (emissiveFactorLocation >= 0) {
        glUniform3f(emissiveFactorLocation, (float)material.emissiveFactor[0],
            (float)material.emissiveFactor[1],
            (float)material.emissiveFactor[2]);
      }
      
      //// Emissive Texture ////
      if (emissiveTextureLocation >= 0) {
        auto textureObject = 0u;
        if (material.emissiveTexture.index >= 0) {
          const auto &texture = model.textures[material.emissiveTexture.index];
          if (texture.source >= 0) {
            textureObject = textureObjects[texture.source];
          }
        }
        glActiveTexture(GL_TEXTURE2);
        glBindTexture(GL_TEXTURE_2D, textureObject);
        glUniform1i(emissiveTextureLocation, 2);
      }

      //// Occlusions Strenght ////
      if (occlusionStrengthLocation >= 0) {
        glUniform1f(occlusionStrengthLocation, (float)material.occlusionTexture.strength);
      }

      //// Occlusions Texture ////
      if (occlusionTextureLocation >= 0) {
        auto textureObject = whiteTexture;
        if (material.occlusionTexture.index >= 0) {
          const auto &texture = model.textures[material.occlusionTexture.index];
          if (texture.source >= 0) {
            textureObject = textureObjects[texture.source];
          }
        }
        glActiveTexture(GL_TEXTURE3);
        glBindTexture(GL_TEXTURE_2D, textureObject);
        glUniform1i(occlusionTextureLocation, 3);
      }

      //// Normal Map ////
      if(normalMapLocation >= 0) {
        auto textureObject = 0u;
        if (material.normalTexture.index >= 0) {
          const auto &texture = model.textures[material.normalTexture.index];
          if (texture.source >= 0) {
            textureObject = textureObjects[texture.source];
          }
        }
        glActiveTexture(GL_TEXTURE4);
        glBindTexture(GL_TEXTURE_2D, textureObject);
        glUniform1i(normalMapLocation, 4);
      }

    } else {
      // Apply default material
      // Defined here:
      // https://github.com/KhronosGroup/glTF/blob/master/specification/2.0/README.md#reference-material
      // https://github.com/KhronosGroup/glTF/blob/master/specification/2.0/README.md#reference-pbrmetallicroughness3

      //// Base Color Factor ////
      if (baseColorFactorLocation >= 0) {
        glUniform4f(baseColorFactorLocation, 1, 1, 1, 1);
      }

      //// Base Color Texture ////
      if (baseColorTextureLocation >= 0) {
        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_2D, whiteTexture);
        glUniform1i(baseColorTextureLocation, 0);
      }

      //// Metallic Factor ////
      if (metallicFactorLocation >= 0) {
        glUniform1f(metallicFactorLocation, 1.f);
      }

      //// Roughness Factor ////
      if (roughnessFactorLocation >= 0) {
        glUniform1f(roughnessFactorLocation, 1.f);
      }

      //// Metallic Roughness Texture ////
      if (metallicRoughnessTextureLocation >= 0) {
        glActiveTexture(GL_TEXTURE1);
        glBindTexture(GL_TEXTURE_2D, 0);
        glUniform1i(metallicRoughnessTextureLocation, 1);
      }

      //// Emissive Factor ////
      if (emissiveFactorLocation >= 0) {
        glUniform3f(emissiveFactorLocation, 0.f, 0.f, 0.f);
      }
      
      //// Emissive Texture ////
      if (emissiveTextureLocation >= 0) {
        glActiveTexture(GL_TEXTURE2);
        glBindTexture(GL_TEXTURE_2D, 0);
        glUniform1i(emissiveTextureLocation, 2);
      }

      //// Occlusions Strenght ////
      if (occlusionStrengthLocation >= 0) {
        glUniform1f(occlusionStrengthLocation, 0.f);
      }

      //// Occlusions Texture ////
      if (occlusionTextureLocation >= 0) {
        glActiveTexture(GL_TEXTURE3);
        glBindTexture(GL_TEXTURE_2D, 0);
        glUniform1i(occlusionTextureLocation, 3);
      }

      //// Normal Map ////
      if(normalMapLocation >= 0) {
        glActiveTexture(GL_TEXTURE4);
        glBindTexture(GL_TEXTURE_2D, 0);
        glUniform1i(occlusionTextureLocation, 4);
      }
    }

    //// Bool ////
    if(viewNormalMapLocation >= 0) {
      glUniform1i(viewNormalMapLocation, viewNormalMap);
    }
  };

  // Lambda function to draw the scene
  const auto drawScene = [&](const Camera &camera) {
    glViewport(0, 0, m_nWindowWidth, m_nWindowHeight);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    const auto viewMatrix = camera.getViewMatrix();

    //// Lught Direction ////
    if (lightDirectionLocation >= 0) {
      if (lightFromCamera) {
        glUniform3f(lightDirectionLocation, 0, 0, 1);
      } else {
        const auto lightDirectionInViewSpace = glm::normalize(glm::vec3(viewMatrix * glm::vec4(lightDirection, 0.)));
        glUniform3f(lightDirectionLocation, lightDirectionInViewSpace[0], lightDirectionInViewSpace[1], lightDirectionInViewSpace[2]);
      }
    }

    //// Lught Intensity ////
    if (lightIntensityLocation >= 0) {
      glUniform3f(lightIntensityLocation, lightIntensity[0], lightIntensity[1], lightIntensity[2]);
    }

    //// Apply Occlusion ////
    if (applyOcclusionLocation >= 0) {
      glUniform1i(applyOcclusionLocation, applyOcclusion);
    }

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
            glUniformMatrix4fv(modelMatrixLocation,         1, GL_FALSE, glm::value_ptr(modelMatrix)); 
            
            
            const auto &mesh = model.meshes[node.mesh];
            const auto &vaoRange = meshToVertexArrays[node.mesh];



            for (auto primIdx = 0; primIdx < mesh.primitives.size(); primIdx++) {
              const auto vao = vertexArrayObjects[vaoRange.begin + primIdx];
              const auto &primitive = mesh.primitives[primIdx];

              bindMaterial(primitive.material);

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
      ImGui::Text("Application average %.3f ms/frame (%.1f FPS)", 1000.0f / ImGui::GetIO().Framerate, ImGui::GetIO().Framerate);
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

      if (ImGui::CollapsingHeader("Light", ImGuiTreeNodeFlags_DefaultOpen)) {
        static float lightTheta = 0.f;
        static float lightPhi = 0.f;

        if (ImGui::SliderFloat("theta", &lightTheta, 0, glm::pi<float>()) || ImGui::SliderFloat("phi", &lightPhi, 0, 2.f * glm::pi<float>())) {
          const auto sinPhi   = glm::sin(lightPhi);
          const auto cosPhi   = glm::cos(lightPhi);
          const auto sinTheta = glm::sin(lightTheta);
          const auto cosTheta = glm::cos(lightTheta);
          lightDirection      = glm::vec3(sinTheta * cosPhi, cosTheta, sinTheta * sinPhi);
        }

        static glm::vec3 lightColor(1.f, 1.f, 1.f);
        static float lightIntensityFactor = 1.f;

        if (ImGui::ColorEdit3("color", (float *) &lightColor) || ImGui::InputFloat("intensity", &lightIntensityFactor)) {
          lightIntensity = lightColor * lightIntensityFactor;
        }
        ImGui::Checkbox("light from camera", &lightFromCamera);
        ImGui::Checkbox("apply occlusion", &applyOcclusion);
        ImGui::Checkbox("view normal map", &viewNormalMap);
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
