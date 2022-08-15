#include <NBody/Visualization/GLLoader.hpp>
#include <NBody/Visualization/Renderer.hpp>

namespace NBody::Visualization {

void Renderer::draw(
    const VertexArray& va,
    const IndexBuffer& ib,
    const Shader& shader
) {
  shader.bind();
  va.bind();
  ib.bind();
  glDrawElements(
      GL_TRIANGLES,
      ib.getCount(),
      ib.getType(),
      nullptr
  );
}

} // NBody::Visualization