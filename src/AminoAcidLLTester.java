import org.junit.Test;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;

public class AminoAcidLLTester {
    @Test
    // this test was to see if my code is functioning properly by recognizing how the combination of each letter is and
    // if it still would recognized the stop combinations of letters.
    public void aminoSequence1(){
        AminoAcidLL head = AminoAcidLL.createFromRNASequence("GCUACGGAGCUUCGGAGCUAG");
        AminoAcidLL.printLinkedList(head);
        char[] expected = {'A','T','E','L','R','S'};
        char[] actual = head.aminoAcidList();

        assertArrayEquals(expected, actual);
    }
    @Test
    // this test was to see if it would work for a smalller combination and if there was a stop combination inbetween
    // the sequence  which it should not.
    public void aminoSequence2(){
        AminoAcidLL head = AminoAcidLL.createFromRNASequence("UCGUUGUAGAUGGUU");
        AminoAcidLL.printLinkedList(head);
        char[] expected = {'S','L'};
        char[] actual = head.aminoAcidList();

        assertArrayEquals(expected, actual);
    }
    @Test
    // this sequence is to see if it would pass if its in a different expected order and if it fail then it should be working
    // correctly
    public void aminoSequence3(){
        AminoAcidLL head = AminoAcidLL.createFromRNASequence("UCGUUGAUGGUU");
        AminoAcidLL.printLinkedList(head);
        char[] expected = {'M','V','S','L'};
        char[] actual = head.aminoAcidList();

        assertArrayEquals(expected, actual);
    }
    @Test
    // this test was to see if it would work for a single combination
    public void aminoSequence4(){
        AminoAcidLL head = AminoAcidLL.createFromRNASequence("UCG");
        AminoAcidLL.printLinkedList(head);
        char[] expected = {'S'};
        char[] actual = head.aminoAcidList();

        assertArrayEquals(expected, actual);
    }
    @Test
    // this test was to see if it would work if the sequence is the same letters of characters.
    public void aminoSequence5(){
        AminoAcidLL head = AminoAcidLL.createFromRNASequence("UGUUGUUGUUGUUGUUGUUGUUGUUGUUGUUGU");
        AminoAcidLL.printLinkedList(head);
        char[] expected = {'C'};
        char[] actual = head.aminoAcidList();

        assertArrayEquals(expected, actual);
    }
    @Test
    // this is to Test the  aminoAcidCompare method to see if it is working properly
    public void aminoAcidCompare1(){
        AminoAcidLL a = AminoAcidLL.createFromRNASequence("GAGGAGACCACCUGCGACUAG");
        AminoAcidLL b = AminoAcidLL.createFromRNASequence("GGUGGUGAGGAGGAGACCACCUAG");
        a = AminoAcidLL.sort(a);
        b = AminoAcidLL.sort(b);
        assertEquals(7, a.aminoAcidCompare(b));
    }
    @Test
    // this test is to seee if the rna sequence method is working
    public void rnaSequence1(){
        AminoAcidLL head = AminoAcidLL.createFromRNASequence("CUAGCACCGCUU");
        AminoAcidLL sorted = AminoAcidLL.sort(head);
        int[] expected = {1,2,1};
        int[] actual = sorted.aminoAcidCounts();

        assertArrayEquals(expected, actual);
    }
    @Test
    // Tests to see isSorted method with a sorted list of sequenced char.
    public void isSorted(){
        AminoAcidLL test = AminoAcidLL.createFromRNASequence("GCGGCAGCCGCUUGGGAC");
        assertEquals(false, test.isSorted());
    }
    @Test
    // Tests to see if it would pass it all amino acids sequences were in that A section.
    public void isSorted2(){
        AminoAcidLL test = AminoAcidLL.createFromRNASequence("GCAGCCGCGGCU");
        assertEquals(true, test.isSorted());
    }
    @Test
    // Tests isSorted method to see if all Stop sequence would still allow it to pass and is in order.
    public void isSorted3(){
        AminoAcidLL test = AminoAcidLL.createFromRNASequence("UGDAUAGUAA");
        assertEquals(true, test.isSorted());
    }










}
