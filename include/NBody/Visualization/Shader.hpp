#ifndef CPP_PROJECT_INCLUDE_NBODY_VISUALIZATION_SHADER_HPP_
#define CPP_PROJECT_INCLUDE_NBODY_VISUALIZATION_SHADER_HPP_

#include <NBody/Visualization/GLLoader.hpp>
#include <string>
#include <unordered_map>
#include <vector>

namespace NBody::Visualization {

/// \brief The abstraction of attribute location in OpenGL.
struct AttribLocation {
  /// \brief The location of the attribute.
  GLint location;
  /// \brief The name of the attribute.
  std::basic_string<GLchar> name;
};

/// \brief The abstraction of shader and program in OpenGL.
class Shader {
 private:
  GLuint program_;
  std::unordered_map<std::string, GLint> locationMap_;

 public:
  /// \brief Constructs a shader with vertex shader source, fragment shader
  /// source and attribute locations. Meanwhile, binds it to current OpenGL
  /// context
  /// \note Runtime exception will be thrown if the construction process of the
  /// shader failed. For example, the OpenGL program cannot compile the shader
  /// source due to syntax error.
  /// \param vShaderSource
  /// \param fShaderSource
  /// \param attribLocations
  Shader(
      const std::basic_string<GLchar>& vShaderSource,
      const std::basic_string<GLchar>& fShaderSource,
      const std::vector<AttribLocation>& attribLocations
  );

  /// \brief Destructs a shader. Automatically deletes the OpenGL program.
  ~Shader();

  /// \brief Binds the shader.
  void bind() const;

  /// \brief Unbinds the shader on the current OpenGL context.
  static void unbind();

  /// \brief Sets the uniform1f variable.
  /// \param name The name of the variable  in shader source
  /// \param f The float value
  void setUniform1f(const std::string& name, GLfloat f);

  /// \brief Sets the uniform4f variable.
  /// \param name The name of the variable in shader source
  /// \param f0 The first float value
  /// \param f1 The second float value
  /// \param f2 The third float value
  /// \param f3 The fourth float value
  void setUniform4f(
      const std::string& name,
      GLfloat f0,
      GLfloat f1,
      GLfloat f2,
      GLfloat f3
  );

  /// \brief Sets the uniformMatrix4fv variable.
  /// \param name The name of the variable in shader source
  /// \param value The pointer to matrix data will be loaded into the uniform
  /// variable
  /// \param transpose Whether to transpose the matrix before loading into the
  /// uniform variable
  void setUniformMatrix4fv(
      const std::string& name,
      const GLfloat* value,
      bool transpose = false
  );

 private:
  GLint getUniformLocation(const std::string& name);
  static GLuint compileShader(
      GLuint type,
      const std::basic_string<GLchar>& source
  );
};

} // NBody::Visualization

#endif // CPP_PROJECT_INCLUDE_NBODY_VISUALIZATION_SHADER_HPP_
