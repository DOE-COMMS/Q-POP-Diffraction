from __future__ import annotations

from dataclasses import dataclass
from html import escape
from pathlib import Path
import re
import shutil


ROOT = Path(__file__).resolve().parent
BUILD_DIR = ROOT / "_build" / "html"
STYLE_FILE = ROOT / "style.css"


@dataclass
class Page:
    source: Path
    title: str
    body: str
    output_name: str


HEADING_RE = re.compile(r"^(#{1,6})\s+(.*)$")
LIST_RE = re.compile(r"^-\s+(.*)$")
LINK_RE = re.compile(r"\[([^\]]+)\]\(([^)]+)\)")
INLINE_CODE_RE = re.compile(r"`([^`]+)`")


def slug_from_name(path: Path) -> str:
    return f"{path.stem}.html"


def apply_inline_markup(text: str) -> str:
    escaped = escape(text)
    escaped = INLINE_CODE_RE.sub(lambda match: f"<code>{escape(match.group(1))}</code>", escaped)

    def replace_link(match: re.Match[str]) -> str:
        label = escape(match.group(1))
        target = match.group(2)
        if target.endswith(".md"):
            target = f"{Path(target).stem}.html"
        return f'<a href="{escape(target, quote=True)}">{label}</a>'

    return LINK_RE.sub(replace_link, escaped)


def render_markdown(text: str) -> tuple[str, str]:
    lines = text.splitlines()
    blocks: list[str] = []
    paragraph: list[str] = []
    list_items: list[str] = []
    code_lines: list[str] = []
    code_lang: str | None = None
    title = "Q-POP-Diffraction Docs"

    def flush_paragraph() -> None:
        nonlocal paragraph
        if paragraph:
            content = " ".join(line.strip() for line in paragraph)
            blocks.append(f"<p>{apply_inline_markup(content)}</p>")
            paragraph = []

    def flush_list() -> None:
        nonlocal list_items
        if list_items:
            items = "".join(f"<li>{apply_inline_markup(item)}</li>" for item in list_items)
            blocks.append(f"<ul>{items}</ul>")
            list_items = []

    def flush_code() -> None:
        nonlocal code_lines, code_lang
        if code_lang is None:
            return
        code = "\n".join(code_lines)
        if code_lang == "mermaid":
            blocks.append(f'<pre class="mermaid">{escape(code)}</pre>')
        else:
            class_attr = f' class="language-{escape(code_lang, quote=True)}"' if code_lang else ""
            blocks.append(f"<pre><code{class_attr}>{escape(code)}</code></pre>")
        code_lines = []
        code_lang = None

    for line in lines:
        if line.startswith("```"):
            flush_paragraph()
            flush_list()
            if code_lang is None:
                code_lang = line[3:].strip()
                code_lines = []
            else:
                flush_code()
            continue

        if code_lang is not None:
            code_lines.append(line)
            continue

        if not line.strip():
            flush_paragraph()
            flush_list()
            continue

        heading_match = HEADING_RE.match(line)
        if heading_match:
            flush_paragraph()
            flush_list()
            level = len(heading_match.group(1))
            heading_text = heading_match.group(2).strip()
            if level == 1:
                title = heading_text
            blocks.append(f"<h{level}>{apply_inline_markup(heading_text)}</h{level}>")
            continue

        list_match = LIST_RE.match(line)
        if list_match:
            flush_paragraph()
            list_items.append(list_match.group(1).strip())
            continue

        paragraph.append(line)

    flush_paragraph()
    flush_list()
    flush_code()
    return title, "\n".join(blocks)


# Navigation order; pages not listed follow alphabetically.
ORDER = ["index", "installation", "inputs", "method", "multi-gpu", "outputs", "testing", "diffraction-flow"]


def load_pages() -> list[Page]:
    pages: list[Page] = []
    rank = {name: i for i, name in enumerate(ORDER)}
    for source in sorted(ROOT.glob("*.md"), key=lambda path: (rank.get(path.stem, len(ORDER)), path.stem)):
        title, body = render_markdown(source.read_text(encoding="utf-8"))
        pages.append(Page(source=source, title=title, body=body, output_name=slug_from_name(source)))
    return pages


def render_template(page: Page, pages: list[Page]) -> str:
    nav_links = []
    for item in pages:
        class_name = "current" if item.output_name == page.output_name else ""
        nav_links.append(
            f'<a class="nav-link {class_name}" href="{escape(item.output_name, quote=True)}">{escape(item.title)}</a>'
        )
    navigation = "".join(nav_links)
    return f"""<!DOCTYPE html>
<html lang=\"en\">
<head>
  <meta charset=\"utf-8\">
  <meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">
  <title>{escape(page.title)}</title>
  <link rel=\"stylesheet\" href=\"style.css\">
  <script type=\"module\">
    import mermaid from 'https://cdn.jsdelivr.net/npm/mermaid@11/dist/mermaid.esm.min.mjs';
    mermaid.initialize({{ startOnLoad: true, theme: 'neutral', securityLevel: 'loose' }});
  </script>
</head>
<body>
  <div class=\"site-shell\">
    <aside class=\"sidebar\">
      <div class=\"sidebar-inner\">
        <p class=\"eyebrow\">Q-POP-Diffraction</p>
        <h1 class=\"site-title\">Implementation Notes</h1>
        <nav>{navigation}</nav>
      </div>
    </aside>
    <main class=\"content\">
      {page.body}
    </main>
  </div>
</body>
</html>
"""


def build_site() -> None:
    if BUILD_DIR.exists():
        shutil.rmtree(BUILD_DIR)
    BUILD_DIR.mkdir(parents=True, exist_ok=True)

    shutil.copy2(STYLE_FILE, BUILD_DIR / "style.css")
    pages = load_pages()

    for page in pages:
        output_path = BUILD_DIR / page.output_name
        output_path.write_text(render_template(page, pages), encoding="utf-8")


if __name__ == "__main__":
    build_site()