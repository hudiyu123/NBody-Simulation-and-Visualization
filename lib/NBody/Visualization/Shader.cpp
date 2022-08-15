#include <NBody/Visualization/Shader.hpp>
#include <stdexcept>

namespace NBody::Visualization {

Shader::Shader(
    const std::basic_string<GLchar>& vShaderSource,
    const std::basic_string<GLchar>& fShaderSource,
    const std::vector<AttribLocation>& attribLocations
) {
  GLuint vShader = compileShader(GL_VERTEX_SHADER, vShaderSource);
  if (!vShader) {
    throw std::runtime_error("Cannot compile vertex shader!");
  }
  GLuint fShader = compileShader(GL_FRAGMENT_SHADER, fShaderSource);
  if (!fShader) {
    glDeleteShader(vShader);
    throw std::runtime_error("Cannot compile fragment shader!");
  }
  GLuint program = glCreateProgram();
  if (!program) {
    throw std::runtime_error("Cannot create program!");
  }
  GLint status = GL_FALSE;
  glAttachShader(program, vShader);
  glAttachShader(program, fShader);
  for (const auto& attribLocation: attribLocations) {
    glBindAttribLocation(
        program,
        attribLocation.location,
        attribLocation.name.c_str()
    );
  }
  glLinkProgram(program);
  glGetProgramiv(program, GL_LINK_STATUS, &status);
  glDeleteShader(vShader);
  glDeleteShader(fShader);
  if (status != GL_TRUE) {
    glDeleteProgram(program);
    throw std::runtime_error("Cannot link program!");
  }
  program_ = program;
}

Shader::~Shader() {
  glDeleteProgram(program_);
}

void Shader::bind() const {
  glUseProgram(program_);
}

void Shader::unbind() {
  glUseProgram(0);
}

void Shader::setUniform1f(const std::string& name, GLfloat f) {
  auto location = getUniformLocation(name);
  glUniform1f(location, f);
}

void Shader::setUniform4f(
    const std::string& name,
    GLfloat f0,
    GLfloat f1,
    GLfloat f2,
    GLfloat f3
) {
  auto location = getUniformLocation(name);
  glUniform4f(location, f0, f1, f2, f3);
}
void Shader::setUniformMatrix4fv(const std::string& name,
                                 const GLfloat* value,
                                 bool transpose) {
  auto location = getUniformLocation(name);
  glUniformMatrix4fv(location, 1, transpose, value);
}

GLint Shader::getUniformLocation(const std::string& name) {
  auto element = locationMap_.find(name);
  if (element != locationMap_.cend()) {
    return element->second;
  } else {
    GLint location = glGetUniformLocation(program_, name.c_str());
    locationMap_[name] = location;
    return location;
  }
}

GLuint Shader::compileShader(
    GLuint type,
    const std::basic_string<GLchar>& source
) {
  GLuint shader = glCreateShader(type);
  if (!shader) {
    return 0;
  }
  const GLchar* cp = &source[0];
  auto len = static_cast<GLint>(source.size());
  glShaderSource(shader, 1, &cp, &len);
  glCompileShader(shader);
  GLint status = GL_FALSE;
  glGetShaderiv(shader, GL_COMPILE_STATUS, &status);
  if (status != GL_TRUE) {
    glDeleteShader(shader);
    return 0;
  }
  return shader;
}

} // NBody::Visualization