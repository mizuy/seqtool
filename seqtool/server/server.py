from bottle import route, run, template, request, response

DEBUG = True



@route('/')
def index():
    return template('''
<h1>seqtool</h1>
<ul>
    <li><a href="seqview">seqview</a></li>
    <li><a href="tssview">tssview</a></li>
    <li><a href="geneview">geneview</a></li>
</ul>
''')


@route('/seqview')
@route('/seqview', method='POST')
def seqview():
    if request.method =='POST':
        seqv = request.forms.get('seqv')
    else:
        seqv = request.cookies.get('seqv', '')
    response.set_cookie('seqv', seqv)
    return template('''
<h1>seqview</h1>
<form action="/seqview" method="post">
  <textarea name="seqv" cols="40", rows="10">{{seqv}}</textarea>
  <input type="submit" value="Submit" />
</form>
<pre>
{{seqv}}
</pre>
''', seqv=seqv)

@route('/tssview')
def tssview():
    return template('''
<h1>tssview</h1>
''')

@route('/geneview')
def geneview():
    return template('''
<h1>geneview</h1>
hoge
''')


def main():
    run(host='localhost', port=8080, debug=DEBUG, reloader=DEBUG)
