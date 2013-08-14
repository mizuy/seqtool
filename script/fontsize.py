import seqtool.svg as svg

def main():
    f = svg.SvgItemsVStack()
    gen = svg.SvgItemGenerator(1,1)

    fonts = ["Andale Mono", "Courier", "Courier New", "Menlo", "Monaco", "Osaka-mono"]

    for font in fonts:
        for fontsize in xrange(5, 30):
            t = svg.SvgItems()
            text = "Sample Text by %s, size=%s" % (font, fontsize)
            fw = fontsize * 0.6
            fh = fontsize * 1.1
            w = fw*len(text)
            h = fh
            t.add(gen.rect(x=0, y=0, width=w, height=h, stroke='red', style='fill:none;', **{'stroke-width': 0.1}))
            t.add(gen.text(text, 0, 0, fontsize=fontsize, font=font))
            f.add(gen.rect(x=0, y=0, width=1, height=5))
            f.add(t)

    print f.svg(f.rect.width, f.rect.height)

if __name__=='__main__':
    main()
