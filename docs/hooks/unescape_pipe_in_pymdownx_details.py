def on_page_markdown(markdown, page,**kwargs):
    out = []
    for line in markdown.splitlines(True):
        # Only touch PyMdown Blocks fence lines
        sline = line.strip()
        if sline.startswith("///"):
            line = line.replace(r"\|", "|")
        elif sline.startswith("![") and sline.endswith(")") and (']{' in sline):
            # Find the indices of the brackets
            line = line.replace(']{','](')
            line = line.replace(')','}')
            line = line.replace('}(','){')
        elif sline.startswith("\clearboth"):
            line = '<br style="clear: both;" />\n'
        out.append(line)
    return "".join(out)
