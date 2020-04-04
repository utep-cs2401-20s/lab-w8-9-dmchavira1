class AminoAcidLL {
    char aminoAcid;
    String[] codons;
    int[] counts;
    AminoAcidLL next;

    AminoAcidLL() {
    }

    /********************************************************************************************/
   AminoAcidLL(String inCodon) {
        this.aminoAcid = AminoAcidResources.getAminoAcidFromCodon(inCodon);
        this.codons = AminoAcidResources.getCodonListForAminoAcid(this.aminoAcid);
        this.counts = new int[codons.length];
        incrementCount(inCodon);
        this.next = null;
    }

    /********************************************************************************************/
    /* Recursive method that increments the count for a specific codon:
     * If it should be at this node, increments it and stops,
     * if not passes the task to the next node.
     * If there is no next node, add a new node to the list that would contain the codon.
     */
    private void addCodon(String inCodon) {
        if (aminoAcid == AminoAcidResources.getAminoAcidFromCodon(inCodon)) {
            this.incrementCount(inCodon);
        }
       else if (next != null) {
            next.addCodon(inCodon);
        }
        else {
            this.next = new AminoAcidLL(inCodon);
        }
    }

    /********************************************************************************************/
    /* Helper method that increments the count of the codon usage */
    private void incrementCount(String inCodon) {
        for (int i = 0; i < this.codons.length; i++) {
            if (this.codons[i].equals(inCodon)) {
                this.counts[i]++;
            }
        }
    }

    /********************************************************************************************/
    /* Shortcut to find the total number of instances of this amino acid */
    private int totalCount() {
        int sum = 0;
        for (int i = 0; i < counts.length; i++) {
            sum += counts[i];
        }
        return sum;
    }

    /********************************************************************************************/
    /* helper method for finding the list difference on two matching nodes
     *  must be matching, but this is not tracked */
    private int totalDiff(AminoAcidLL inList) {
        return Math.abs(totalCount() - inList.totalCount());
    }

    /********************************************************************************************/
    /* helper method for finding the list difference on two matching nodes
     *  must be matching, but this is not tracked */
    private int codonDiff(AminoAcidLL inList) {
        int diff = 0;
        for (int i = 0; i < codons.length; i++) {
            diff += Math.abs(counts[i] - inList.counts[i]);
        }
        return diff;
    }

    /********************************************************************************************/
    /* Recursive method that finds the differences in **Amino Acid** counts.
     * the list *must* be sorted to use this method */
    public int aminoAcidCompare(AminoAcidLL inList) {
        if(inList.next == null && this.next == null){
            if(this.aminoAcid == inList.aminoAcid)
                return this.totalDiff(inList);
            return inList.totalCount() + this.totalCount();
        }

        if(this.aminoAcid == inList.aminoAcid){
            if(this.next == null) {
                return this.totalDiff(inList) + sum(inList.next);
            }
            if(inList.next == null)
                return this.totalDiff(inList) + sum(this.next);

            return this.totalDiff(inList) + this.next.aminoAcidCompare(inList.next);
        }

       if(this.next != null && this.aminoAcid < inList.aminoAcid){
            return this.totalCount() + this.next.aminoAcidCompare(inList);
        }

        return inList.totalCount() + this.aminoAcidCompare(inList.next);
    }
    /********************************************************************************************/
    /* Same as above, but counts the codon usage differences
     * Must be sorted. */
    public int codonCompare(AminoAcidLL inList) {
        if (inList.next == null) {
            if (this.next == null){
                if (this.aminoAcid == inList.aminoAcid)
                    return this.codonDiff(inList);

            return inList.totalCount() + this.totalCount();
        }
    }

        if(this.aminoAcid == inList.aminoAcid){
           if(this.next == null) {
                return this.codonDiff(inList) + sum(inList.next);
            }
            if(inList.next == null)
                return this.codonDiff(inList) + sum(this.next);
            return this.codonDiff(inList) + this.next.codonCompare(inList.next);
        }

        if(this.aminoAcid < inList.aminoAcid) {
            if (this.next != null) {
                return this.totalCount() + this.next.codonCompare(inList);
            }
        }
        return inList.totalCount() + this.codonCompare(inList.next);
    }


    /********************************************************************************************/
    /* Recursively returns the total list of amino acids in the order that they are in in the linked list. */
    public char[] aminoAcidList() {
        if (next == null) {
            return new char[]{aminoAcid};
        }
        char[] temp = next.aminoAcidList();
       char[] ret = new char[temp.length+1];
        ret[0] = aminoAcid;
        for(int i = 0; i < temp.length; i++){
            ret[i+1] = temp[i];
        }
        return ret;
    }
    /********************************************************************************************/
    /* Recursively returns the total counts of amino acids in the order that they are in in the linked list. */
    public int[] aminoAcidCounts() {
        if(next == null){
            return new int[]{this.totalCount()};
        }
       int[] temp = next.aminoAcidCounts();
        int[] ret = new int[temp.length+1];
          ret[0] = totalCount();
        for(int i = 0; i < temp.length; i++){
            ret[i+1] = temp[i];
        }
        return ret;
    }
    public static void printAminoAcidCounts(int[] array){
        for(int i = 0; i < array.length; i++){
            System.out.print(array[i]);
        }
        System.out.println();
    }

    /********************************************************************************************/
    /* recursively determines if a linked list is sorted or not */
    public boolean isSorted() {
        if(next == null){
            return true;
        }
        if(aminoAcid > next.aminoAcid)
            return false;
        return next.isSorted();
    }
    /********************************************************************************************/
    /* Static method for generating a linked list from an RNA sequence */
    public static AminoAcidLL createFromRNASequence(String inSequence) {
       AminoAcidLL head = null;

        for (int i = 0; i < inSequence.length(); i++) {
            if(AminoAcidResources.getAminoAcidFromCodon(inSequence.substring(0,3)) == '*')
                break;
            if (i == 0) {
                head = new AminoAcidLL(inSequence.substring(0, 3));
                inSequence = inSequence.substring(3);
            }
            else {
                if(AminoAcidResources.getAminoAcidFromCodon(inSequence.substring(0,3)) == '*')
                    break;
                head.addCodon(inSequence.substring(0, 3));
                inSequence = inSequence.substring(3);
                i = 1;
            }
        }
        return head;
    }

    /********************************************************************************************/
    /* sorts a list by amino acid character*/
    public static AminoAcidLL sort(AminoAcidLL inList) {
        AminoAcidLL beforeCurrent = inList;
        AminoAcidLL curNode = inList.next;
        AminoAcidLL next;
        AminoAcidLL position;

        while(curNode != null){
            next = curNode.next;
            position = findInsertion(curNode.aminoAcid, inList);
            if(position == beforeCurrent){
                beforeCurrent = curNode;
            }
            else {
                beforeCurrent.next = null;
                if(position == null){
                    AminoAcidLL temp = inList;
                    inList = curNode;
                    curNode.next = temp;
                    beforeCurrent.next = next;
                }
                else{
                    position.next = curNode;
                    curNode.next = beforeCurrent;
                    beforeCurrent.next = next;
                }
            }
            curNode = next;
        }
        return inList;
    }

    /* helper method that searches the list from the head toward the current node to find the insertion position */
    private static AminoAcidLL findInsertion(char aminoAcid, AminoAcidLL inList){
        AminoAcidLL curNodeA = null;
        AminoAcidLL curNodeB = inList;

        while(aminoAcid > curNodeB.aminoAcid && curNodeB != null){
            curNodeA = curNodeB;
            curNodeB = curNodeB.next;
        }
        return curNodeA;
    }

    /********************************************************************************************/
    /* helper method which allows me to see if the code is functioning well*/
    public static void printLinkedList(AminoAcidLL head){
        AminoAcidLL temp = head;

        System.out.println(" The order of linked list is:");
        while(temp != null){
            System.out.print(temp.aminoAcid + " ");
            temp = temp.next;
        }
        System.out.println();
    }
    //helper method
    public int sum(AminoAcidLL node){
        AminoAcidLL temp = node;
        int sum = 0;
        while(temp  != null){
            sum += temp.totalCount();
            temp = temp.next;
        }

        return sum;
    }
}