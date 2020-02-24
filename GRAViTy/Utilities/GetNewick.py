def GetNewick(node, newick, parentdist, leaf_names):
	if node.is_leaf():
		return "%s:%f%s" % (leaf_names[node.id], parentdist - node.dist, newick)
	else:
		if len(newick) > 0:
			newick = "):%f%s" % (parentdist - node.dist, newick)
		else:
			newick = ");"
		newick = GetNewick(node.get_left(), newick, node.dist, leaf_names)
		newick = GetNewick(node.get_right(), ",%s" % (newick), node.dist, leaf_names)
		newick = "(%s" % (newick)
		return newick

