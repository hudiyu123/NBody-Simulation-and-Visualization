#ifndef CPP_PROJECT_INCLUDE_NBODY_VISUALIZATION_GL_LOADER_HPP_
#define CPP_PROJECT_INCLUDE_NBODY_VISUALIZATION_GL_LOADER_HPP_

#ifdef __APPLE__
#include <GL/glew.h>
#else
#include <GL/gl.h>

namespace NBody::Visualization {

#define PAPAYA_GL_LIST \
GLE(void,      AttachShader,            GLuint program, GLuint shader) \
GLE(void,      BindAttribLocation,      GLuint program, GLuint index, const GLchar* name) \
GLE(void,      BindBuffer,              GLenum target, GLuint buffer) \
GLE(void,      BindVertexArray,         GLuint array) \
GLE(void,      BufferData,              GLenum target, GLsizeiptr size, const GLvoid* data, GLenum usage) \
GLE(void,      BufferSubData,           GLenum target, GLintptr offset, GLsizeiptr size, const GLvoid* data) \
GLE(void,      CompileShader,           GLuint shader) \
GLE(GLuint,    CreateProgram,           void) \
GLE(GLuint,    CreateShader,            GLenum type) \
GLE(void,      DeleteBuffers,           GLsizei n, const GLuint* buffers) \
GLE(void,      DeleteShader,            GLuint shader) \
GLE(void,      DeleteProgram,           GLuint program) \
GLE(void,      DeleteVertexArrays,      GLsizei n, const GLuint* arrays) \
GLE(void,      DrawBuffers,             GLsizei n, const GLenum* bufs) \
GLE(void,      Enable,                  GLenum cap) \
GLE(void,      EnableVertexAttribArray, GLuint index) \
GLE(void,      GenVertexArrays,         GLsizei n, GLuint* arrays) \
GLE(void,      GenBuffers,              GLsizei n, GLuint* buffers) \
GLE(void,      GetProgramiv,            GLuint program, GLenum pname, GLint* params) \
GLE(GLint,     GetAttribLocation,       GLuint program, const GLchar* name) \
GLE(void,      GetShaderiv,             GLuint shader, GLenum pname, GLint* params) \
GLE(GLint,     GetUniformLocation,      GLuint program, const GLchar* name) \
GLE(void,      LinkProgram,             GLuint program) \
GLE(void,      ShaderSource,            GLuint shader, GLsizei count, const GLchar* const* string, const GLint* length) \
GLE(void,      Uniform1f,               GLint location, GLfloat v0) \
GLE(void,      Uniform4f,               GLint location, GLfloat v0, GLfloat v1, GLfloat v2, GLfloat v3) \
GLE(void,      UniformMatrix4fv,        GLint location, GLsizei count, GLboolean transpose, const GLfloat* value) \
GLE(void,      UseProgram,              GLuint program) \
GLE(void,      VertexAttribPointer,     GLuint index, GLint size, GLenum type, GLboolean normalized, GLsizei stride, const GLvoid*  pointer) \

#define GLE(ret, name, ...) \
typedef ret name##proc(__VA_ARGS__); \
extern name##proc* gl##name;
PAPAYA_GL_LIST
#undef GLE

/// \brief Initializes the GL Loader.
/// \return true if the initialization succeeded, false otherwise
bool glLoaderInit();

} // NBody::Visualization
#endif

#endif // CPP_PROJECT_INCLUDE_NBODY_VISUALIZATION_GL_LOADER_HPP_
