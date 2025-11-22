/*
 * Copyright 2024-2025, Seqera Labs
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package nextflow.script.control;

import org.codehaus.groovy.ast.ASTNode;
import org.codehaus.groovy.control.SourceUnit;
import org.codehaus.groovy.control.messages.WarningMessage;
import org.codehaus.groovy.syntax.CSTNode;

/**
 * A warning that should only be reported when paranoid warnings
 * are enabled.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
public class ParanoidWarning extends WarningMessage implements RelatedInformationAware {

    private String otherMessage;

    private ASTNode otherNode;

    public ParanoidWarning(int importance, String message, CSTNode context, SourceUnit owner) {
        super(importance, message, context, owner);
    }

    public void setRelatedInformation(String otherMessage, ASTNode otherNode) {
        this.otherMessage = otherMessage;
        this.otherNode = otherNode;
    }

    @Override
    public String getOtherMessage() {
        return otherMessage;
    }

    @Override
    public ASTNode getOtherNode() {
        return otherNode;
    }
}
