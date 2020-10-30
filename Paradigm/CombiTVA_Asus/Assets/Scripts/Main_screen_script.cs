using System.Collections;
using System.Collections.Generic;
using UnityEngine;
//L added
using UnityEngine.UI;
using System.Linq; // added
using System.Threading;

// Written for Nexus 10 tablet
// Lulu 190724 Unity version 2019
// 190725 version without coroutine StartParadigm
// 190905 Specialise for Nexus

public class Main_screen_script : MonoBehaviour
{
    // On screen objects
    public Button startbutton, confirmbutton;
    public Text titletext, starttext, responsetextUI;
    public Dropdown[] dropdownsUI;
    public InputField[] inputsUI;
    public Image redcross;
    public Image[] maskimagesUI;
    public Text[] stimtextUI;
    public Texture backsp, enter;

    // Script objects
    int nruns,currun,curtrial,pos;
    float expdur,actualexpdur,currCountdownValue;
    List<int> triallist;
    bool iskBoardOpen,confirmresp,gettingresp,isdisplaying;//isandroid, 
    Coroutine wf;

    // Experiment display objects
    Sprite[] maskimagesbmp;
    int curstate,stimdur,correctresponses, allresponses;
    int[] stimdurarray = {1, 2, 3, 5, 9, 12};
    int[] maskind = { 0, 1, 2, 3, 4, 5, 6, 7 };
    int[] letterind, trialarray;
    int[][] stimpos_unilat = { 
        new int[] {0,1}, new int[] { 0, 2 }, new int[] { 1, 2 }, 
        new int[] { 3, 4 }, new int[] { 3, 5 }, new int[] { 4, 5 } 
    };
    int[][][] stimpos_bilat = {
        new int[][] {new int[] {0,3}, new int[] { 0, 4 }, new int[] { 0,5 }},
        new int[][] {new int[] {1,3}, new int[] { 1, 4 }, new int[] { 1,5 }},
        new int[][] {new int[] {2,3}, new int[] { 2, 4 }, new int[] { 2,5 }}
    };
    string curtrialstate,ips, savedata,targetstring,distractorstring,filename,instructionstxt, pc;
    string[] letters = { "A", "B", "D", "E", "F", "G", "H", "J", "K", "L", "M", "N", "O", "P", "R", "S", "T", "V", "X", "Z" };
    
    // -------------- Start experiment --------------
    void Start()
    {
        // Timing
        QualitySettings.vSyncCount = 1;
        float Rate = 60.0f;
        Application.targetFrameRate = (int)Rate*10;
        Cursor.visible = true;
        
        // Display
        int xRes = 1920; int yRes = 1080; Screen.SetResolution(xRes, yRes, true, (int)Rate); // resolution x, y, fullscreen, refresh rate
        int viewDis = 600; float screenSize = 609.6f;
        float posRdeg = 7.5f;
        float posradius = Mathf.Tan(posRdeg / 180 * Mathf.PI) * viewDis / (screenSize / Mathf.Sqrt(xRes * xRes + yRes * yRes));
        float[] posangles = { Mathf.PI / 3f, 0f, -Mathf.PI / 3f, -2 * Mathf.PI / 3f, Mathf.PI, 2 * Mathf.PI / 3f };
        for (int i = 0; i < 6; i++)
        {
            Vector3 stimposition = new Vector3(posradius * Mathf.Cos(posangles[i]), posradius * Mathf.Sin(posangles[i]), 0);
            maskimagesUI[i].GetComponent<RectTransform>().localPosition = stimposition;
            maskimagesUI[i].gameObject.SetActive(false);
            stimtextUI[i].GetComponent<RectTransform>().localPosition = stimposition;
            stimtextUI[i].gameObject.SetActive(false);
        }
        maskimagesbmp = Resources.LoadAll<Sprite>("Masks");
        redcross.gameObject.SetActive(false);
        responsetextUI.gameObject.SetActive(false);
        startbutton.gameObject.SetActive(false);
        starttext.gameObject.SetActive(false);

        // Exp settings
        //isandroid = true;
        instructionstxt = starttext.text;
        curstate = 0; curtrial = 0; currun = 0; correctresponses = 0; allresponses = 0;
        letterind = Enumerable.Range(0, 20).ToArray();
        curtrialstate = "preptrial";
        iskBoardOpen = false; confirmresp = false; gettingresp = false;
    }

    private void Update()
    {
        switch (curstate)
        {
            case 0: // on hold
                if (Input.GetKeyDown(KeyCode.Return)) { TaskOnClick(); }
                if (Input.GetKeyDown(KeyCode.Escape)) { Application.Quit(); }
                break;
            case 1: //start paradigm
                StopCoroutine(StartCountdown());
                goto case 2;
            case 2: // task ongoing until end is reached
                if (curtrial == trialarray.Length)
                {
                    isdisplaying = false;
                    System.IO.File.AppendAllText(filename, "#End of run " + (currun + 1).ToString() + " Percentage correct: " + (100 * correctresponses / allresponses).ToString("F0") + "% " + Time.time.ToString("F3") + "s\n");
                    curstate = 3; break;//goto case 3;
                }
                switch (curtrialstate)
                {
                    case "preptrial":
                        stimdur = PrepareTrial(trialarray[curtrial], currun);
                        curtrialstate = "showfixation";responsetextUI.text = "";
                        break;
                    case "showfixation":
                        redcross.gameObject.SetActive(true);
                        savedata = "\t"+Time.time.ToString("F4");
                        wf = StartCoroutine(WFrames(60));
                        isdisplaying = true; curtrialstate = "showstim";
                        break;
                    case "showstim":
                        if (isdisplaying == false)
                        {
                            StopCoroutine(wf);
                            savedata = savedata+"\t"+Time.time.ToString("F4");
                            redcross.gameObject.SetActive(true);
                            for (int i = 0; i < stimtextUI.Length; i++){stimtextUI[i].gameObject.SetActive(true);}
                            wf = StartCoroutine(WFrames(stimdur));
                            isdisplaying = true;  curtrialstate = "showmask";
                        }
                        break;
                    case "showmask":
                        if (isdisplaying == false) //false)//
                        {
                            StopCoroutine(wf);
                            savedata = savedata+"\t"+Time.time.ToString("F4");
                            actualexpdur = expdur;
                            for (int i = 0; i < stimtextUI.Length; i++)
                            {
                                stimtextUI[i].gameObject.SetActive(false);
                                maskimagesUI[i].gameObject.SetActive(true);
                            }
                            wf = StartCoroutine(WFrames(30));
                            isdisplaying = true; curtrialstate = "report";
                        }
                        break;
                    case "report":
                        if (isdisplaying == false) //false)//
                        {
                            StopCoroutine(wf); savedata = savedata+"\t"+Time.time.ToString("F4")+"\tkp:";
                            for (int i = 0; i < stimtextUI.Length; i++){maskimagesUI[i].gameObject.SetActive(false);}
                            redcross.gameObject.SetActive(false);
                            responsetextUI.gameObject.SetActive(true);
                            gettingresp = true; isdisplaying = true;
                        }
                        if (gettingresp) { 
                            string ipstmp = Input.inputString;
                            foreach (char c in ipstmp)
                            {
                                savedata = savedata + c.ToString() + "-" + Time.time.ToString("F4") + ";";
                                if (c == '\b' && responsetextUI.text.Length != 0) // has backspace/delete been pressed?
                                {
                                    responsetextUI.text = responsetextUI.text.Substring(0, responsetextUI.text.Length - 1);
                                }
                                else if ((c == '\n') || (c == '\r')) // enter/return
                                {
                                    confirmresp = true; gettingresp = false;
                                }
                                else if (letters.Contains(c.ToString().ToUpper()) == true && responsetextUI.text.Contains(c.ToString().ToUpper()) == false && responsetextUI.text.Length < 6)
                                {
                                    responsetextUI.text += c.ToString().ToUpper();
                                }
                            }
                        }
                        if (confirmresp)
                        {
                            ips = responsetextUI.text;
                            responsetextUI.gameObject.SetActive(false);
                            if (ips.Length == 0){
                                ips="-";
                            }
                            else {
                                allresponses += ips.Length;
                                foreach (char r in ips) { if (targetstring.Contains(r.ToString())) { correctresponses += 1; } }
                            }
                            savedata = (1 + trialarray[curtrial]).ToString() + "\t" + (1000 * actualexpdur).ToString("F1") + "\t" + targetstring + "\t" + distractorstring + "\t" + ips + savedata + "\n";
                            System.IO.File.AppendAllText(filename, savedata);
                            curtrial++;
                            confirmresp = false; gettingresp = false; 
                            curtrialstate = "preptrial";
                        }
                        break;
                }
                break;
            case 3:
                // feedback
                if (isdisplaying == false) {
                    starttext.gameObject.SetActive(true);
                    if (allresponses > 0) {
                        pc = (100 * correctresponses / allresponses).ToString("F0") + "%";
                        if ((100 * correctresponses / allresponses)>90){ pc = pc + "\n\nTry to report more letters and feel free to guess more"; }
                        else if ((100 * correctresponses / allresponses)<80){  pc = pc + "\n\nTry to make fewer guesses and report only the letters you are relatively sure of"; }
                    }
                    else { pc = "0"; }
                    if (dropdownsUI[3].value==1){
                        starttext.text = "End of run "+(currun+1).ToString()+" of 9\nPercentage correct: " + pc;}
                    else {starttext.text = "End of the practice run\nPercentage correct: " + pc;}
                    isdisplaying = true;
                }
                if (Input.GetKeyDown(KeyCode.Return) || Input.GetMouseButtonDown(0))
                {
                    if (currun < nruns-1)
                    {
                        currun += 1; curtrial = 0; correctresponses = 0; allresponses = 0;
                        curstate = 0;
                        StartCoroutine(StartCountdown());
                    }
                    else
                    {
                        goto case 4;
                    }
                }
                break;
            case 4:
            responsetextUI.gameObject.SetActive(true);
            if (nruns>1){
                starttext.text = "End\n\nThank you for your participation!";
                curstate = 4;
                if (Input.GetKeyDown(KeyCode.Return)){Application.Quit();}
            } else {
                for (int i = 0; i<dropdownsUI.Length;i++){ dropdownsUI[i].gameObject.SetActive(true); }
                for (int i = 0; i<inputsUI.Length;i++){ inputsUI[i].gameObject.SetActive(true); }
                confirmbutton.gameObject.SetActive(true);titletext.gameObject.SetActive(true);
                starttext.text = instructionstxt; starttext.gameObject.SetActive(true);
                dropdownsUI[3].value=1;
                Start();
            }
            break;
        }
    }


    // -------------- Initialisation functions (start button) --------------
    public void ConfirmSettings(){
        for (int i = 0; i<dropdownsUI.Length;i++){ dropdownsUI[i].gameObject.SetActive(false); }
        for (int i = 0; i<inputsUI.Length;i++){ inputsUI[i].gameObject.SetActive(false); }
        confirmbutton.gameObject.SetActive(false);
        titletext.gameObject.SetActive(false);
        startbutton.gameObject.SetActive(true);
        starttext.gameObject.SetActive(true);
        
        string partSex = dropdownsUI[0].options[dropdownsUI[0].value].text;
        string partHand = dropdownsUI[1].options[dropdownsUI[1].value].text;
        int partSession = dropdownsUI[2].value;
        int selectedmode = dropdownsUI[3].value;
        int selectedkeyboard = dropdownsUI[4].value;
        string partID = inputsUI[0].text; if (partID.Length < 1) { partID = "0"; }
        string partAge = inputsUI[1].text; string partEdu = inputsUI[2].text;
        
        triallist = Enumerable.Range(0, 24).ToList();
        Text startbuttontext = startbutton.gameObject.GetComponentInChildren(typeof(Text)) as Text;
        if (selectedmode > 0) {
            for (int i = 0; i < 6; i++) { triallist.Add(i); triallist.Add(i); }
            startbuttontext.text = "Press enter to start experiment";
        } else{
            startbuttontext.text = "Press enter to start practice";
        }
        nruns = (selectedmode * 8) + 1; trialarray = triallist.ToArray(); Shuffle(trialarray);
        filename = "P" + partID + "_s" + (partSession+1).ToString() + "D-" + selectedmode.ToString() + ".txt";
        if (Application.isMobilePlatform){filename = Application.persistentDataPath + "/" + filename;}
        System.IO.File.WriteAllText(filename, "#"+System.DateTime.Now.ToString("yyMMdd_HHmm")+ "\tCombiTVA - " + startbuttontext.text.Substring(21,startbuttontext.text.Length-21)+"\n#Participant ID/age/sex/handedness/education/session: "+partID+"/"+partAge+"/"+ partSex+"/"+partHand+"/"+ partEdu+ "/" +(partSession+1).ToString()+"\n");
        System.IO.File.AppendAllText(filename, "#Desktop version - Stim positions:\t");
        for (int i = 0; i< 6; i++){System.IO.File.AppendAllText(filename, maskimagesUI[i].GetComponent<RectTransform>().localPosition.ToString());}
        System.IO.File.AppendAllText(filename, "\n#Resolution, refresh rate:"+Screen.currentResolution.ToString());

        System.IO.File.AppendAllText(filename, "\n"+ (nruns * trialarray.Length).ToString() + "\n");
    }
    
    public void TaskOnClick() {
        startbutton.gameObject.SetActive(false);
        Cursor.visible = false;
        StartCoroutine(StartCountdown());
    }

    // -------------- Timing stuff -------------- 
    public IEnumerator StartCountdown(int countdownValue = 3)
    {
        starttext.gameObject.SetActive(true);
        currCountdownValue = countdownValue;
        while (currCountdownValue > 0)
        {
            starttext.text = currCountdownValue.ToString();
            yield return new WaitForSeconds(1.0f);
            currCountdownValue--;
        }
        starttext.gameObject.SetActive(false);
        System.IO.File.AppendAllText(filename, "\n#Start run "+ (currun+1).ToString() + " "+ Time.time.ToString("F3")+ "s\n");
        curstate = 1; 
    }

    public IEnumerator WFrames(int frameCount) //static
    {
        expdur = 0;
        while (frameCount > 0)
        {
            frameCount--;
            yield return new WaitForEndOfFrame();
            expdur += Time.smoothDeltaTime;
        }
        isdisplaying = false; 
    }

    // -------------- Paradigm support functions --------------
    public int PrepareTrial(int trialtype, int run)
    {
        distractorstring = ""; targetstring = "";
        Shuffle(maskind); Shuffle(letterind);
        for (int i = 0; i < 6; i++)
        {
            maskimagesUI[i].sprite = maskimagesbmp[maskind[i]];
            // first fill all letters, then remove or change colors
            stimtextUI[i].text = letters[letterind[i]];
            stimtextUI[i].color = Color.red;
            if (trialtype > 5) // not whole report
            {
                if (trialtype < 12)
                { // whole report 2T unilateral
                    pos = System.Array.IndexOf(stimpos_unilat[trialtype - 6], i);
                    if (pos < 0) { stimtextUI[i].text = ""; targetstring = targetstring + "0"; } //distractorstring = distractorstring + "0"; 
                    else { targetstring = targetstring + letters[letterind[i]]; }
                    distractorstring = distractorstring + "0";
                }
                else if (trialtype < 15)
                { //whole report 2T bilateral
                    pos = System.Array.IndexOf(stimpos_bilat[run % 3][trialtype - 12], i);
                    if (pos < 0) { stimtextUI[i].text = ""; targetstring = targetstring + "0"; }
                    else { targetstring = targetstring + letters[letterind[i]]; }
                    distractorstring = distractorstring + "0";
                }
                else if (trialtype < 21)
                { // partial report unilateral
                    pos = System.Array.IndexOf(stimpos_unilat[trialtype - 15], i);
                    if (pos < 0)
                    {
                        stimtextUI[i].color = Color.blue; distractorstring = distractorstring + letters[letterind[i]];
                        targetstring = targetstring + "0";
                    }
                    else { targetstring = targetstring + letters[letterind[i]]; distractorstring = distractorstring + "0"; }
                }
                else
                { // partial report bilateral
                    pos = System.Array.IndexOf(stimpos_bilat[(run + 1) % 3][trialtype - 21], i);
                    if (pos < 0)
                    {
                        stimtextUI[i].color = Color.blue; distractorstring = distractorstring + letters[letterind[i]];
                        targetstring = targetstring + "0";
                    }
                    else { targetstring = targetstring + letters[letterind[i]]; distractorstring = distractorstring + "0"; }
                }
                stimdur = stimdurarray[3];
            }
            else
            {
                targetstring = targetstring + letters[letterind[i]];
                distractorstring = distractorstring + "0";
                stimdur = stimdurarray[trialarray[curtrial]]; 
            }
        }
        return stimdur;
    }

    void Shuffle(int[] array)
    {
        System.Random random = new System.Random();
        int n = array.Count();
        while (n > 1)
        {
            n--;
            int i = random.Next(n + 1);
            int temp = array[i];
            array[i] = array[n];
            array[n] = temp;
        }
    }
}
