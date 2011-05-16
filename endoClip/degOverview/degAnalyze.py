
import bioLibCG

def getExpressionByType(fN, eCutoff):

        type_expr = {}
        type_count = {}

        f = open(fN, 'r')
        for line in f:
                ls = line.strip().split('\t')
                expr = int(ls[3])
                type = ls[7]
                
                if not expr >= int(eCutoff): continue
                type_expr[type] = type_expr.get(type, 0) + expr
                type_count[type] = type_count.get(type, 0) + 1

        for type, expr in type_expr.iteritems():
                print type, expr, type_count[type]

if __name__ == "__main__":
        import sys
        bioLibCG.submitArgs(globals()[sys.argv[1]], sys.argv[1:])
