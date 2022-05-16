#!/usr/bin/env python3

# Simple library to read a Newick tree

import sys


# Object to use to read a string character by character
class StringReader:

    def __init__(self, string):
        # Only need to keep a reference to the string and our position
        self.string = string
        self.pos = 0

    # Have we reached the end of the string
    def isExhausted(self):
        return self.pos == len(self.string)

    # Return the character at the current position
    def getCurrentCharacter(self):
        return self.string[self.pos]

    # Go to the next position
    def advance(self):
        self.pos += 1

    # Move until a character from `characters` is found
    def readUntil(self, characters):
        startPos = self.pos
        while not self.isExhausted():
            if self.getCurrentCharacter() in characters:
                return self.string[startPos:self.pos]
            self.advance()
        return self.string[startPos:]


# Function to build a tree from a Newick string
# The tree is returned as a dictionary that has two keys:
#  - label: the name of the node
#  - children: the list of the children, themselves represented as dictionaries
# A leaf has no children (empty list)
def buildTree(tree_string: str):
    treeReader = StringReader(tree_string)

    def _rec_buildTree():
        if treeReader.getCurrentCharacter() == '(':
            # Opening bracket -> a sub-tree of at least 1 child, separated by commas
            treeReader.advance()
            children = []
            # Read the first child
            firstChild = _rec_buildTree()
            children.append(firstChild)
            # Read more children as long as separated by a comma
            while treeReader.getCurrentCharacter() == ',':
                treeReader.advance()
                nextChild = _rec_buildTree()
                children.append(nextChild)
            # We should be hitting the closing bracket
            assert treeReader.getCurrentCharacter() == ')', treeReader.string[:treeReader.pos]
            treeReader.advance()
            # Here comes the label. NB: NHX tags are included in the label
            label = treeReader.readUntil('),')
            node = {'label': label, 'children': children}
            return node
        else:
            # A leaf
            label = treeReader.readUntil('),')
            node = {'label': label}
            return node
    return _rec_buildTree()


if __name__ == "__main__":
    s = sys.stdin.readline()
    t = buildTree(s)
    print(t)
